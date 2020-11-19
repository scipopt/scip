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
/*  along with SCIP; see the file COPYING. If not visit scipopt.org.         */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   cons_nonlinear.c
 * @ingroup DEFPLUGINS_CONS
 * @brief  constraint handler for nonlinear constraints specified by algebraic expressions
 * @author Ksenia Bestuzheva
 * @author Benjamin Mueller
 * @author Felipe Serrano
 * @author Stefan Vigerske
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

#include "scip/cons_nonlinear.h"
#include "scip/expr_var.h"
#include "scip/expr_sum.h"
#include "scip/cons_linear.h"
#include "scip/cons_and.h"
#include "scip/cons_bounddisjunction.h"
#include "scip/heur_subnlp.h"
#include "scip/heur_trysol.h"
#include "scip/debug.h"
#include "scip/struct_nlhdlr.h"

/* fundamental constraint handler properties */
#define CONSHDLR_NAME          "nonlinear"
#define CONSHDLR_DESC          "handler for nonlinear constraints specified by algebraic expressions"
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

/* properties of the nonlinear constraint handler statistics table */
#define TABLE_NAME_NONLINEAR           "nonlinear"
#define TABLE_DESC_NONLINEAR           "nonlinear constraint handler statistics"
#define TABLE_POSITION_NONLINEAR       12500                  /**< the position of the statistics table */
#define TABLE_EARLIEST_STAGE_NONLINEAR SCIP_STAGE_TRANSFORMED /**< output of the statistics table is only printed from this stage onwards */

#define VERTEXPOLY_MAXPERTURBATION      1e-3 /**< maximum perturbation */
#define VERTEXPOLY_USEDUALSIMPLEX       TRUE /**< use dual or primal simplex algorithm? */
#define VERTEXPOLY_RANDNUMINITSEED  20181029 /**< seed for random number generator, which is used to move points away from the boundary */
#define VERTEXPOLY_ADJUSTFACETFACTOR     1e1 /**< adjust resulting facets in checkRikun() up to a violation of this value times lpfeastol */

#define BRANCH_RANDNUMINITSEED      20191229 /**< seed for random number generator, which is used to select from several similar good branching candidates */

#define BILIN_MAXNAUXEXPRS                10 /**< maximal number of auxiliary expressions per bilinear term */

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

/** enforcement data of an expression */
typedef struct
{
   SCIP_NLHDLR*          nlhdlr;             /**< nonlinear handler */
   SCIP_NLHDLREXPRDATA*  nlhdlrexprdata;     /**< data of nonlinear handler */
   SCIP_CONSNONLINEAR_EXPRENFO_METHOD nlhdlrparticipation; /**< methods where nonlinear handler participates */
   SCIP_Bool             issepainit;         /**< was the initsepa callback of nlhdlr called */
   SCIP_Real             auxvalue;           /**< auxiliary value of expression w.r.t. currently enforced solution */
   SCIP_Bool             sepabelowusesactivity;/**< whether sepabelow uses activity of some expression */
   SCIP_Bool             sepaaboveusesactivity;/**< whether sepaabove uses activity of some expression */
} EXPRENFO;

/** data stored by constraint handler in an expression that belongs to a nonlinear constraint */
struct SCIP_Expr_OwnerData
{
   SCIP_CONSHDLR*        conshdlr;           /** nonlinear constraint handler */

   /* locks */
   int                   nlockspos;          /**< positive locks counter */
   int                   nlocksneg;          /**< negative locks counter */

   /* propagation (in addition to activity that is stored in expr) */
   SCIP_INTERVAL         propbounds;         /**< bounds to propagate in reverse propagation */
   unsigned int          propboundstag;      /**< tag to indicate whether propbounds are valid for the current propagation rounds */
   SCIP_Bool             inpropqueue;        /**< whether expression is queued for propagation */

   /* enforcement of expr == auxvar (or expr <= auxvar, or expr >= auxvar) */
   EXPRENFO**            enfos;              /**< enforcements */
   int                   nenfos;             /**< number of enforcements, or -1 if not initialized */
   unsigned int          lastenforced;       /**< last enforcement round where expression was enforced successfully */
   unsigned int          nactivityusesprop;  /**< number of nonlinear handlers whose activity computation (or domain propagation) depends on the activity of the expression */
   unsigned int          nactivityusessepa;  /**< number of nonlinear handlers whose separation (estimate or enfo) depends on the activity of the expression */
   unsigned int          nauxvaruses;        /**< number of nonlinear handlers whose separation uses an auxvar in the expression */
   SCIP_VAR*             auxvar;             /**< auxiliary variable used for outer approximation cuts */

   /* branching */
   SCIP_Real             violscoresum;       /**< sum of violation scores for branching stored for this expression */
   SCIP_Real             violscoremax;       /**< max of violation scores for branching stored for this expression */
   int                   nviolscores;        /**< number of violation scores stored for this expression */
   unsigned int          violscoretag;       /**< tag to decide whether a violation score of an expression needs to be initialized */

   /* additional data for variable expressions (TODO move into sub-struct?) */
   SCIP_CONS**           conss;              /**< constraints in which this variable appears */
   int                   nconss;             /**< current number of constraints in conss */
   int                   consssize;          /**< length of conss array */
   SCIP_Bool             consssorted;        /**< is the array of constraints sorted */

   int                   filterpos;          /**< position of eventdata in SCIP's event filter, -1 if not catching events */
};

/** constraint data for nonlinear constraints */
struct SCIP_ConsData
{
   /* data that defines the constraint: expression and sides */
   SCIP_EXPR*            expr;               /**< expression that represents this constraint */
   SCIP_Real             lhs;                /**< left-hand side */
   SCIP_Real             rhs;                /**< right-hand side */

   /* variables */
   SCIP_EXPR**           varexprs;           /**< array containing all variable expressions */
   int                   nvarexprs;          /**< total number of variable expressions */
   SCIP_Bool             catchedevents;      /**< do we catch events on variables? */

   /* constraint violation */
   SCIP_Real             lhsviol;            /**< violation of left-hand side by current solution */
   SCIP_Real             rhsviol;            /**< violation of right-hand side by current solution */
   SCIP_Real             gradnorm;           /**< norm of gradient of constraint function in current solution (if evaluated) */
   unsigned int          gradnormsoltag;     /**< tag of solution used that gradnorm corresponds to */

   /* status flags */
   unsigned int          ispropagated:1;     /**< did we propagate the current bounds already? */
   unsigned int          issimplified:1;     /**< did we simplify the expression tree already? */

   /* locks */
   int                   nlockspos;          /**< number of positive locks */
   int                   nlocksneg;          /**< number of negative locks */

   /* repair infeasible solutions */
   SCIP_VAR*             linvardecr;         /**< variable that may be decreased without making any other constraint infeasible, or NULL if none */
   SCIP_VAR*             linvarincr;         /**< variable that may be increased without making any other constraint infeasible, or NULL if none */
   SCIP_Real             linvardecrcoef;     /**< linear coefficient of linvardecr */
   SCIP_Real             linvarincrcoef;     /**< linear coefficient of linvarincr */

   /* miscellaneous */
   SCIP_EXPRCURV         curv;               /**< curvature of the root expression w.r.t. the original variables */
   SCIP_NLROW*           nlrow;              /**< a nonlinear row representation of this constraint */
   int                   consindex;          /**< an index of the constraint that is unique among all expr-constraints in this SCIP instance and is constant */
};

/** constraint upgrade method */
typedef struct
{
   SCIP_DECL_NONLINCONSUPGD((*consupgd));    /**< method to call for upgrading nonlinear constraint */
   int                   priority;           /**< priority of upgrading method */
   SCIP_Bool             active;             /**< is upgrading enabled */
} CONSUPGRADE;

/** constraint handler data */
struct SCIP_ConshdlrData
{
   /* nonlinear handler */
   SCIP_NLHDLR**         nlhdlrs;            /**< nonlinear handlers */
   int                   nnlhdlrs;           /**< number of nonlinear handlers */
   int                   nlhdlrssize;        /**< size of nlhdlrs array */
   SCIP_Bool             indetect;           /**< whether we are currently in detectNlhdlr */
   SCIP_Bool             registerusesactivitysepabelow; /**< a flag that is used only during \ref @detectNlhdlr() */
   SCIP_Bool             registerusesactivitysepaabove; /**< a flag that is used only during \ref @detectNlhdlr() */

   /* constraint upgrades */
   CONSUPGRADE**         consupgrades;       /**< constraint upgrade methods for specializing nonlinear constraints */
   int                   consupgradessize;   /**< size of consupgrades array */
   int                   nconsupgrades;      /**< number of constraint upgrade methods */

   /* other plugins */
   SCIP_EVENTHDLR*       eventhdlr;          /**< handler for variable bound change events */
   SCIP_HEUR*            subnlpheur;         /**< a pointer to the subnlp heuristic, if available */
   SCIP_HEUR*            trysolheur;         /**< a pointer to the trysol heuristic, if available */

   /* tags and counters */
   int                   auxvarid;           /**< unique id for the next auxiliary variable */
   unsigned int          curboundstag;       /**< tag indicating current variable bounds */
   unsigned int          lastboundrelax;     /**< tag when bounds where most recently relaxed */
   unsigned int          lastvaractivitymethodchange; /**< tag when method used to evaluate activity of variables changed last */
   unsigned int          enforound;          /**< total number of enforcement calls, including current one */
   int                   lastconsindex;      /**< last used consindex, plus one */

   /* activity intervals and domain propagation */
   SCIP_DECL_EXPR_INTEVALVAR((*intevalvar)); /**< method currently used for activity calculation of variable expressions */
   SCIP_Bool             globalbounds;       /**< whether global variable bounds should be used for activity calculation */
   SCIP_QUEUE*           reversepropqueue;   /**< expression queue to be used in reverse propagation, filled by SCIPtightenExprIntervalNonlinear */
   SCIP_Bool             forceboundtightening; /**< whether bound change passed to SCIPtightenExprIntervalNonlinear should be forced */
   unsigned int          curpropboundstag;   /**< tag indicating current propagation rounds, to match with expr->propboundstag */

   /* parameters */
   int                   maxproprounds;      /**< limit on number of propagation rounds for a set of constraints within one round of SCIP propagation */
   SCIP_Bool             propauxvars;        /**< whether to check bounds of all auxiliary variable to seed reverse propagation */
   char                  varboundrelax;      /**< strategy on how to relax variable bounds during bound tightening */
   SCIP_Real             varboundrelaxamount; /**< by how much to relax variable bounds during bound tightening */
   SCIP_Real             conssiderelaxamount; /**< by how much to relax constraint sides during bound tightening */
   SCIP_Real             vp_maxperturb;      /**< maximal relative perturbation of reference point */
   SCIP_Real             vp_adjfacetthreshold; /**< adjust computed facet up to a violation of this value times lpfeastol */
   SCIP_Bool             vp_dualsimplex;     /**< whether to use dual simplex instead of primal simplex for facet computing LP */
   SCIP_Bool             reformbinprods;     /**< whether to reformulate products of binary variables during presolving */
   SCIP_Bool             reformbinprodsand;  /**< whether to use the AND constraint handler for reformulating binary products */
   int                   reformbinprodsfac;  /**< minimum number of terms to reformulate bilinear binary products by factorizing variables (<= 1: disabled) */
   SCIP_Bool             forbidmultaggrnlvar; /**< whether to forbid multiaggregation of variables that appear in a nonlinear term of a constraint */
   SCIP_Bool             tightenlpfeastol;   /**< whether to tighten LP feasibility tolerance during enforcement, if it seems useful */
   SCIP_Bool             propinenforce;      /**< whether to (re)run propagation in enforcement */
   SCIP_Real             weakcutthreshold;   /**< threshold for when to regard a cut from an estimator as weak */
   SCIP_Real             strongcutmaxcoef;   /**< "strong" cuts will be scaled to have their maximal coef in [1/strongcutmaxcoef,strongcutmaxcoef] */
   SCIP_Bool             strongcutefficacy;  /**< consider efficacy requirement when deciding whether a cut is "strong" */
   SCIP_Bool             forcestrongcut;     /**< whether to force "strong" cuts in enforcement */
   SCIP_Real             enfoauxviolfactor;  /**< an expression will be enforced if the "auxiliary" violation is at least enfoauxviolfactor times the "original" violation */
   SCIP_Real             weakcutminviolfactor; /**< retry with weak cuts for constraints with violation at least this factor of maximal violated constraints */
   char                  violscale;          /**< method how to scale violations to make them comparable (not used for feasibility check) */
   char                  checkvarlocks;      /**< whether variables contained in a single constraint should be forced to be at their lower or upper bounds ('d'isable, change 't'ype, add 'b'ound disjunction) */
   int                   branchauxmindepth;  /**< from which depth on to allow branching on auxiliary variables */
   SCIP_Bool             branchexternal;     /**< whether to use external branching candidates for branching */
   SCIP_Real             branchhighviolfactor; /**< consider a constraint highly violated if its violation is >= this factor * maximal violation among all constraints */
   SCIP_Real             branchhighscorefactor; /**< consider a variable branching score high if its branching score >= this factor * maximal branching score among all variables */
   SCIP_Real             branchviolweight;   /**< weight by how much to consider the violation assigned to a variable for its branching score */
   SCIP_Real             branchdualweight;   /**< weight by how much to consider the dual values of rows that contain a variable for its branching score */
   SCIP_Real             branchpscostweight; /**< weight by how much to consider the pseudo cost of a variable for its branching score */
   SCIP_Real             branchdomainweight; /**< weight by how much to consider the domain width in branching score */
   SCIP_Real             branchvartypeweight;/**< weight by how much to consider variable type in branching score */
   char                  branchscoreagg;     /**< how to aggregate several branching scores given for the same expression ('a'verage, 'm'aximum, or 's'um) */
   char                  branchviolsplit;    /**< method used to split violation in expression onto variables ('e'venly, 'm'idness of solution, 'd'omain width, 'l'ogarithmic domain width) */
   SCIP_Real             branchpscostreliable; /**< minimum pseudo-cost update count required to consider pseudo-costs reliable */

   /* statistics */
   SCIP_Longint          nweaksepa;          /**< number of times we used "weak" cuts for enforcement */
   SCIP_Longint          ntightenlp;         /**< number of times we requested solving the LP with a smaller feasibility tolerance when enforcing */
   SCIP_Longint          ndesperatetightenlp; /**< number of times we requested solving the LP with a smaller feasibility tolerance when enforcing because we didn't know anything better */
   SCIP_Longint          ndesperatebranch;   /**< number of times we branched on some variable because normal enforcement was not successful */
   SCIP_Longint          ndesperatecutoff;   /**< number of times we cut off a node in enforcement because no branching candidate could be found */
   SCIP_Longint          nforcelp;           /**< number of times we forced solving the LP when enforcing a pseudo solution */
   SCIP_CLOCK*           canonicalizetime;   /**< time spend for canonicalization */
   SCIP_Longint          ncanonicalizecalls; /**< number of times we called canonicalization */

   /* facets of envelops of vertex-polyhedral functions */
   SCIP_RANDNUMGEN*      vp_randnumgen;      /**< random number generator used to perturb reference point */
   SCIP_LPI*             vp_lp[SCIP_MAXVERTEXPOLYDIM+1];  /**< LPs used to compute facets for functions of different dimension */

   /* hashing of bilinear terms */
   SCIP_HASHTABLE*       bilinhashtable;     /**< hash table for bilinear terms */
   SCIP_CONSNONLINEAR_BILINTERM* bilinterms; /**< bilinear terms */
   int                   nbilinterms;        /**< total number of bilinear terms */
   int                   bilintermssize;     /**< size of bilinterms array */
   int                   bilinmaxnauxexprs;  /**< maximal number of auxiliary expressions per bilinear term */

   /* branching */
   SCIP_RANDNUMGEN*      branchrandnumgen;   /**< random number generated used in branching variable selection */
   char                  branchpscostupdatestrategy; /**< value of parameter branching/lpgainnormalize */

   /* misc */
   SCIP_Bool             checkedvarlocks;    /**< whether variables contained in a single constraint have been already considered */
   SCIP_HASHMAP*         var2expr;           /**< hashmap to map SCIP variables to variable-expressions */
};

/** branching candidate with various scores */
typedef struct
{
   SCIP_EXPR*            expr;               /**< expression that holds branching candidate */
   SCIP_Real             auxviol;            /**< aux-violation score of candidate */
   SCIP_Real             domain;             /**< domain score of candidate */
   SCIP_Real             dual;               /**< dual score of candidate */
   SCIP_Real             pscost;             /**< pseudo-cost score of candidate */
   SCIP_Real             vartype;            /**< variable type score of candidate */
   SCIP_Real             weighted;           /**< weighted sum of other scores, see scoreBranchingCandidates() */
} BRANCHCAND;

/*
 * Local methods
 */

/** frees auxiliary variables of expression, if any */
static
SCIP_RETCODE freeAuxVar(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EXPR*            expr                /**< expression which auxvar to free, if any */
   )
{
   SCIP_EXPR_OWNERDATA* mydata;

   assert(scip != NULL);
   assert(expr != NULL);

   mydata = SCIPexprGetOwnerData(expr);
   assert(mydata != NULL);

   if( mydata->auxvar == NULL )
      return SCIP_OKAY;

   SCIPdebugMsg(scip, "remove auxiliary variable <%s> for expression %p\n", SCIPvarGetName(mydata->auxvar), (void*)expr);

   /* remove variable locks
    * as this is a relaxation-only variable, no other plugin should use it for deducing any type of reductions or cutting planes
    */
   SCIP_CALL( SCIPaddVarLocks(scip, mydata->auxvar, -1, -1) );

   /* release auxiliary variable */
   SCIP_CALL( SCIPreleaseVar(scip, &mydata->auxvar) );
   assert(mydata->auxvar == NULL);

   return SCIP_OKAY;
}

/** frees data used for enforcement of expression, that is, nonlinear handlers
 *
 * can also clear indicators whether expr needs enforcement methods, that is,
 * free an associated auxiliary variable and reset the activityusage counts
 */
static
SCIP_RETCODE freeEnfoData(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EXPR*            expr,               /**< expression whose enforcement data will be released */
   SCIP_Bool             freeauxvar          /**< whether aux var should be released and activity usage counts be reset */
   )
{
   SCIP_EXPR_OWNERDATA* mydata;
   int e;

   mydata = SCIPexprGetOwnerData(expr);
   assert(mydata != NULL);

   if( freeauxvar )
   {
      /* free auxiliary variable */
      SCIP_CALL( freeAuxVar(scip, expr) );
      assert(mydata->auxvar == NULL);

      /* reset count on activity and auxvar usage */
      mydata->nactivityusesprop = 0;
      mydata->nactivityusessepa = 0;
      mydata->nauxvaruses = 0;
   }
#if !1  // FIXME
   /* free data stored by nonlinear handlers */
   for( e = 0; e < mydata->nenfos; ++e )
   {
      SCIP_NLHDLR* nlhdlr;

      assert(mydata->enfos[e] != NULL);

      nlhdlr = mydata->enfos[e]->nlhdlr;
      assert(nlhdlr != NULL);

      if( mydata->enfos[e]->issepainit )
      {
         /* call the separation deinitialization callback of the nonlinear handler */
         SCIP_CALL( SCIPcallNlhdlrExitSepaNonlinear(scip, nlhdlr, expr, mydata->enfos[e]->nlhdlrexprdata) );
         mydata->enfos[e]->issepainit = FALSE;
      }

      /* free nlhdlr exprdata, if there is any and there is a method to free this data */
      if( mydata->enfos[e]->nlhdlrexprdata != NULL && nlhdlr->freeexprdata != NULL )
      {
         SCIP_CALL( nlhdlr->freeexprdata(scip, nlhdlr, expr, &mydata->enfos[e]->nlhdlrexprdata) );
         assert(mydata->enfos[e]->nlhdlrexprdata == NULL);
      }

      /* free enfo data */
      SCIPfreeBlockMemory(scip, &mydata->enfos[e]);
   }
#endif
   /* free array with enfo data */
   SCIPfreeBlockMemoryArrayNull(scip, &mydata->enfos, mydata->nenfos);

   /* we need to look at this expression in detect again */
   mydata->nenfos = -1;

   return SCIP_OKAY;
}

/** callback that frees data that this conshdlr stored in an expression */
static
SCIP_DECL_EXPR_OWNERDATAFREE(exprownerdataFree)
{
   assert(scip != NULL);
   assert(expr != NULL);
   assert(ownerdata != NULL);
   assert(*ownerdata != NULL);

   /* expression should not be locked anymore */
   assert((*ownerdata)->nlockspos == 0);
   assert((*ownerdata)->nlocksneg == 0);

   SCIP_CALL( freeEnfoData(scip, expr, TRUE) );

   /* expression should not be enforced anymore */
   assert((*ownerdata)->nenfos == 0);
   assert((*ownerdata)->auxvar == NULL);

   if( SCIPisExprVar(scip, expr) )
   {
      SCIP_CONSHDLRDATA* conshdlrdata;
      SCIP_VAR* var;

      /* there should be no constraints left that still use this variable */
      assert((*ownerdata)->nconss == 0);
      /* thus, there should also be no variable event catched (via this exprhdlr) */
      assert((*ownerdata)->filterpos == -1);

      SCIPfreeBlockMemoryArrayNull(scip, &(*ownerdata)->conss, (*ownerdata)->consssize);

      /* update var2expr hashmap in conshdlrdata */
      conshdlrdata = SCIPconshdlrGetData((*ownerdata)->conshdlr);
      assert(conshdlrdata != NULL);

      var = SCIPgetVarExprVar(expr);
      assert(var != NULL);

      /* if no variable-expression stored for var hashmap, then the var hasn't been used in any constraint, so do nothing
       * if variable-expression stored for var is different, then also do nothing
       */
      if( SCIPhashmapGetImage(conshdlrdata->var2expr, var) != (void*)expr )
         return SCIP_OKAY;

      /* remove var -> varexpr map from hashmap */
      SCIP_CALL( SCIPhashmapRemove(conshdlrdata->var2expr, var) );
   }

   SCIPfreeBlockMemory(scip, ownerdata);

   return SCIP_OKAY;
}

/** callback that creates data that this conshdlr wants to store in an expression */
static
SCIP_DECL_EXPR_OWNERDATACREATE(exprownerdataCreate)
{
   assert(scip != NULL);
   assert(expr != NULL);
   assert(ownerdata != NULL);

   SCIP_CALL( SCIPallocClearBlockMemory(scip, ownerdata) );
   (*ownerdata)->nenfos = -1;
   (*ownerdata)->conshdlr = (SCIP_CONSHDLR*)ownerdatacreatedata;

   if( SCIPisExprVar(scip, expr) )
   {
      SCIP_CONSHDLRDATA* conshdlrdata;
      SCIP_VAR* var;

      (*ownerdata)->filterpos = -1;

      /* add to var2expr hashmap if not having expr for var yet */

      conshdlrdata = SCIPconshdlrGetData((*ownerdata)->conshdlr);
      assert(conshdlrdata != NULL);

      var = SCIPgetVarExprVar(expr);

      if( !SCIPhashmapExists(conshdlrdata->var2expr, (void*)var) )
      {
         /* store the variable expression in the hashmap */
         SCIP_CALL( SCIPhashmapInsert(conshdlrdata->var2expr, (void*)var, (void*)expr) );
      }
      else
      {
         /* if expr was just created, then it shouldn't already be stored as image of var */
         assert(SCIPhashmapGetImage(conshdlrdata->var2expr, (void*)var) !=  (void*)expr);
      }
   }

   *ownerdatafree = exprownerdataFree;

   return SCIP_OKAY;
}

/** creates a variable expression or retrieves from hashmap in conshdlr data */
static
SCIP_RETCODE createExprVar(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< nonlinear constraint handler */
   SCIP_EXPR**           expr,               /**< pointer where to store expression */
   SCIP_VAR*             var                 /**< variable to be stored */
   )
{
   assert(conshdlr != NULL);
   assert(expr != NULL);
   assert(var != NULL);

   /* get variable expression representing the given variable if there is one already */
   *expr = (SCIP_EXPR*) SCIPhashmapGetImage(SCIPconshdlrGetData(conshdlr)->var2expr, (void*) var);

   if( *expr == NULL )
   {
      /* create a new variable expression; this also captures the expression */
      SCIP_CALL( SCIPcreateExprVar(scip, expr, var, exprownerdataCreate, (SCIP_EXPR_OWNERDATACREATEDATA*)conshdlr) );
      assert(*expr != NULL);
      /* exprownerdataCreate should have added var->expr to var2expr */
      assert(SCIPhashmapGetImage(SCIPconshdlrGetData(conshdlr)->var2expr, (void*)var) == (void*)*expr);
   }
   else
   {
      /* only capture already existing expr to get a consistent uses-count */
      SCIPcaptureExpr(*expr);
   }

   return SCIP_OKAY;
}

#if !1
static
SCIP_DECL_EXPR_MAPVAR(transformVar)
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
#endif

static
SCIP_DECL_EXPR_MAPEXPR(mapexprvar)
{
   SCIP_CONSHDLR* conshdlr = (SCIP_CONSHDLR*)mapexprdata;

   assert(sourcescip != NULL);
   assert(targetscip != NULL);
   assert(sourceexpr != NULL);
   assert(targetexpr != NULL);
   assert(mapexprdata != NULL);

   /* do not provide map if not variable */
   if( !SCIPisExprVar(sourcescip, sourceexpr) )
   {
      *targetexpr = NULL;
      return SCIP_OKAY;
   }

   SCIP_CALL( createExprVar(targetscip, conshdlr, targetexpr, SCIPgetVarExprVar(sourceexpr)) );

   return SCIP_OKAY;
}

/** interval evaluation of variables as used in bound tightening
 *
 * Returns slightly relaxed local variable bounds of a variable as interval.
 * Does not relax beyond integer values, thus does not relax bounds on integer variables at all.
 */
static
SCIP_DECL_EXPR_INTEVALVAR(intEvalVarBoundTightening)
{
   SCIP_INTERVAL interval;
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_Real lb;
   SCIP_Real ub;

   assert(scip != NULL);
   assert(var != NULL);

   conshdlrdata = (SCIP_CONSHDLRDATA*)intevalvardata;
   assert(conshdlrdata != NULL);

   if( conshdlrdata->globalbounds )
   {
      lb = SCIPvarGetLbGlobal(var);
      ub = SCIPvarGetUbGlobal(var);
   }
   else
   {
      lb = SCIPvarGetLbLocal(var);
      ub = SCIPvarGetUbLocal(var);
   }
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

/** processes variable fixing or bound change event */
static
SCIP_DECL_EVENTEXEC(processVarEvent)
{  /*lint --e{715}*/
   SCIP_EVENTTYPE eventtype;
   SCIP_EXPR* expr;
   SCIP_EXPR_OWNERDATA* ownerdata;

   eventtype = SCIPeventGetType(event);
   assert(eventtype & (SCIP_EVENTTYPE_BOUNDCHANGED | SCIP_EVENTTYPE_VARFIXED));

   assert(eventdata != NULL);
   expr = (SCIP_EXPR*) eventdata;
   assert(SCIPisExprVar(scip, expr));

   SCIPdebugMsg(scip, "  exec event %#x for variable <%s> (local [%g,%g], global [%g,%g])\n", eventtype,
         SCIPvarGetName(SCIPeventGetVar(event)),
         SCIPvarGetLbLocal(SCIPeventGetVar(event)), SCIPvarGetUbLocal(SCIPeventGetVar(event)),
         SCIPvarGetLbGlobal(SCIPeventGetVar(event)), SCIPvarGetUbGlobal(SCIPeventGetVar(event)));

   ownerdata = SCIPexprGetOwnerData(expr);
   assert(ownerdata != NULL);
   /* we only catch varevents for variables in constraints, so there should be constraints */
   assert(ownerdata->nconss > 0);
   assert(ownerdata->conss != NULL);

   /* notify constraints that use this variable expression (expr) to repropagate and possibly resimplify
    * - propagation can only find something new if a bound was tightened
    * - simplify can only find something new if a var is fixed (or maybe a bound is tightened)
    *   and we look at global changes (that is, we are not looking at boundchanges in probing)
    */
   if( eventtype & (SCIP_EVENTTYPE_BOUNDTIGHTENED | SCIP_EVENTTYPE_VARFIXED) )
   {
      SCIP_CONSDATA* consdata;
      int c;

      for( c = 0; c < ownerdata->nconss; ++c )
      {
         assert(ownerdata->conss[c] != NULL);
         consdata = SCIPconsGetData(ownerdata->conss[c]);

         /* if boundtightening, then mark constraints to be propagated again
          * TODO we could try be more selective here and only trigger a propagation if a relevant bound has changed,
          *   that is, we don't need to repropagate x + ... <= rhs if only the upper bound of x has been tightened
          *   the locks could help if they were available on a per-constraint base, but they aren't (and it may not be worth it)
          */
         if( eventtype & SCIP_EVENTTYPE_BOUNDTIGHTENED )
         {
            consdata->ispropagated = FALSE;
            SCIPdebugMsg(scip, "  marked <%s> for propagate\n", SCIPconsGetName(ownerdata->conss[c]));
         }

         /* if still in presolve (but not probing), then mark constraints to be unsimplified */
         if( SCIPgetStage(scip) == SCIP_STAGE_PRESOLVING && !SCIPinProbing(scip) )
         {
            consdata->issimplified = FALSE;
            SCIPdebugMsg(scip, "  marked <%s> for simplify\n", SCIPconsGetName(ownerdata->conss[c]));
         }
      }
   }

   /* update curboundstag, lastboundrelax, and expr activity */
   if( eventtype & SCIP_EVENTTYPE_BOUNDCHANGED )
   {
      SCIP_CONSHDLRDATA* conshdlrdata;
      SCIP_INTERVAL activity;

      conshdlrdata = SCIPconshdlrGetData(ownerdata->conshdlr);
      assert(conshdlrdata != NULL);

      /* increase tag on bounds */
      ++conshdlrdata->curboundstag;
      assert(conshdlrdata->curboundstag > 0);

      /* remember also if we relaxed bounds now */
      if( eventtype & SCIP_EVENTTYPE_BOUNDRELAXED )
         conshdlrdata->lastboundrelax = conshdlrdata->curboundstag;

      /* update the activity of the var-expr here immediately
       * (we could call expr->activity = intevalvar(var, consdhlr) directly, but then the exprhdlr statistics are not updated)
       */
      SCIP_CALL( SCIPexprhdlrIntEvalExpr(scip, expr, &activity, conshdlrdata->intevalvar, conshdlrdata) );
#ifdef DEBUG_PROP
      SCIPdebugMsg(scip, "  var-exprhdlr::inteval = [%.20g, %.20g]\n", activity.inf, activity.sup);
#endif
      SCIPexprSetActivity(expr, activity);
      SCIPexprSetActivityTag(expr, conshdlrdata->curboundstag);
   }

   return SCIP_OKAY;
}

/** creates and captures a nonlinear constraint
 *
 * @attention Use copyexpr=FALSE only if expr is already "owned" by conshdlr
 */
static
SCIP_RETCODE createCons(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler */
   SCIP_CONS**           cons,               /**< pointer to hold the created constraint */
   const char*           name,               /**< name of constraint */
   SCIP_EXPR*            expr,               /**< expression of constraint (must not be NULL) */
   SCIP_Real             lhs,                /**< left hand side of constraint */
   SCIP_Real             rhs,                /**< right hand side of constraint */
   SCIP_Bool             copyexpr,           /**< whether to copy the expression or reuse the given expr (capture it) */
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
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONSDATA* consdata;

   assert(conshdlr != NULL);
   assert(expr != NULL);

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

   if( copyexpr )
   {
      /* copy expression, thereby map variables expressions to already existing variables expressions in var2expr map, or augment var2expr map */
      SCIP_CALL( SCIPduplicateExpr(scip, expr, &consdata->expr, mapexprvar, conshdlr, exprownerdataCreate, (SCIP_EXPR_OWNERDATACREATEDATA*)conshdlr) );
   }
   else
   {
      consdata->expr = expr;
      SCIPcaptureExpr(consdata->expr);
   }
   consdata->lhs = lhs;
   consdata->rhs = rhs;
   consdata->consindex = conshdlrdata->lastconsindex++;
   consdata->curv = SCIP_EXPRCURV_UNKNOWN;

   /* create constraint */
   SCIP_CALL( SCIPcreateCons(scip, cons, name, conshdlr, consdata, initial, separate, enforce, check, propagate,
         local, modifiable, dynamic, removable, FALSE) );

   return SCIP_OKAY;
}

/** print statistics for nonlinear handlers */
static
void printNlhdlrStatistics(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< nonlinear constraint handler */
   FILE*                 file                /**< file handle, or NULL for standard out */
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata;
   int i;

   assert(scip != NULL);
   assert(conshdlr != NULL);

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   SCIPinfoMessage(scip, file, "Nlhdlrs            : %10s %10s %10s %10s %10s %10s %10s %10s %10s %10s %10s %10s %10s\n", "Detects", "EnfoCalls", "#IntEval", "PropCalls", "DetectAll", "Separated", "Cutoffs", "DomReds", "BranchScor", "DetectTime", "EnfoTime", "PropTime", "IntEvalTi");

   for( i = 0; i < conshdlrdata->nnlhdlrs; ++i )
   {
      SCIP_NLHDLR* nlhdlr = conshdlrdata->nlhdlrs[i];
      assert(nlhdlr != NULL);

      /* skip disabled nlhdlr */
      if( !nlhdlr->enabled )
         continue;

      SCIPinfoMessage(scip, file, "  %-17s:", nlhdlr->name);
      SCIPinfoMessage(scip, file, " %10lld", nlhdlr->ndetectionslast);
      SCIPinfoMessage(scip, file, " %10lld", nlhdlr->nenfocalls);
      SCIPinfoMessage(scip, file, " %10lld", nlhdlr->nintevalcalls);
      SCIPinfoMessage(scip, file, " %10lld", nlhdlr->npropcalls);
      SCIPinfoMessage(scip, file, " %10lld", nlhdlr->ndetections);
      SCIPinfoMessage(scip, file, " %10lld", nlhdlr->nseparated);
      SCIPinfoMessage(scip, file, " %10lld", nlhdlr->ncutoffs);
      SCIPinfoMessage(scip, file, " %10lld", nlhdlr->ndomreds);
      SCIPinfoMessage(scip, file, " %10lld", nlhdlr->nbranchscores);
      SCIPinfoMessage(scip, file, " %10.2f", SCIPgetClockTime(scip, nlhdlr->detecttime));
      SCIPinfoMessage(scip, file, " %10.2f", SCIPgetClockTime(scip, nlhdlr->enfotime));
      SCIPinfoMessage(scip, file, " %10.2f", SCIPgetClockTime(scip, nlhdlr->proptime));
      SCIPinfoMessage(scip, file, " %10.2f", SCIPgetClockTime(scip, nlhdlr->intevaltime));
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

   SCIPinfoMessage(scip, file, "Enforce            : %10s %10s %10s %10s %10s %10s\n", "WeakSepa", "TightenLP", "DespTghtLP", "DespBranch", "DespCutoff", "ForceLP");
   SCIPinfoMessage(scip, file, "  consexpr%-9s:", "");
   SCIPinfoMessage(scip, file, " %10lld", conshdlrdata->nweaksepa);
   SCIPinfoMessage(scip, file, " %10lld", conshdlrdata->ntightenlp);
   SCIPinfoMessage(scip, file, " %10lld", conshdlrdata->ndesperatetightenlp);
   SCIPinfoMessage(scip, file, " %10lld", conshdlrdata->ndesperatebranch);
   SCIPinfoMessage(scip, file, " %10lld", conshdlrdata->ndesperatecutoff);
   SCIPinfoMessage(scip, file, " %10lld", conshdlrdata->nforcelp);
   SCIPinfoMessage(scip, file, "\n");
   SCIPinfoMessage(scip, file, "Presolve           : %10s\n", "CanonTime");
   SCIPinfoMessage(scip, file, "  consexpr%-9s:", "");
   SCIPinfoMessage(scip, file, " %10.2f", SCIPgetClockTime(scip, conshdlrdata->canonicalizetime));
   SCIPinfoMessage(scip, file, "\n");
}

/*
 * Callback methods of constraint handler
 */

/* TODO: Implement all necessary constraint handler methods. The methods with #if SCIP_DISABLED_CODE ... #else #define ... are optional */

/** copy method for constraint handler plugins (called when SCIP copies plugins) */
#if SCIP_DISABLED_CODE
static
SCIP_DECL_CONSHDLRCOPY(conshdlrCopyNonlinear)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of nonlinear constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define conshdlrCopyNonlinear NULL
#endif

/** destructor of constraint handler to free constraint handler data (called when SCIP is exiting) */
#if SCIP_DISABLED_CODE
static
SCIP_DECL_CONSFREE(consFreeNonlinear)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of nonlinear constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consFreeNonlinear NULL
#endif


/** initialization method of constraint handler (called after problem was transformed) */
#if SCIP_DISABLED_CODE
static
SCIP_DECL_CONSINIT(consInitNonlinear)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of nonlinear constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consInitNonlinear NULL
#endif


/** deinitialization method of constraint handler (called before transformed problem is freed) */
#if SCIP_DISABLED_CODE
static
SCIP_DECL_CONSEXIT(consExitNonlinear)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of nonlinear constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consExitNonlinear NULL
#endif


/** presolving initialization method of constraint handler (called when presolving is about to begin) */
#if SCIP_DISABLED_CODE
static
SCIP_DECL_CONSINITPRE(consInitpreNonlinear)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of nonlinear constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consInitpreNonlinear NULL
#endif


/** presolving deinitialization method of constraint handler (called after presolving has been finished) */
#if SCIP_DISABLED_CODE
static
SCIP_DECL_CONSEXITPRE(consExitpreNonlinear)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of nonlinear constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consExitpreNonlinear NULL
#endif


/** solving process initialization method of constraint handler (called when branch and bound process is about to begin) */
#if SCIP_DISABLED_CODE
static
SCIP_DECL_CONSINITSOL(consInitsolNonlinear)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of nonlinear constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consInitsolNonlinear NULL
#endif


/** solving process deinitialization method of constraint handler (called before branch and bound process data is freed) */
#if SCIP_DISABLED_CODE
static
SCIP_DECL_CONSEXITSOL(consExitsolNonlinear)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of nonlinear constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consExitsolNonlinear NULL
#endif


/** frees specific constraint data */
#if SCIP_DISABLED_CODE
static
SCIP_DECL_CONSDELETE(consDeleteNonlinear)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of nonlinear constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consDeleteNonlinear NULL
#endif


/** transforms constraint data into data belonging to the transformed problem */
#if SCIP_DISABLED_CODE
static
SCIP_DECL_CONSTRANS(consTransNonlinear)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of nonlinear constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consTransNonlinear NULL
#endif


/** LP initialization method of constraint handler (called before the initial LP relaxation at a node is solved) */
#if SCIP_DISABLED_CODE
static
SCIP_DECL_CONSINITLP(consInitlpNonlinear)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of nonlinear constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consInitlpNonlinear NULL
#endif


/** separation method of constraint handler for LP solutions */
#if SCIP_DISABLED_CODE
static
SCIP_DECL_CONSSEPALP(consSepalpNonlinear)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of nonlinear constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consSepalpNonlinear NULL
#endif


/** separation method of constraint handler for arbitrary primal solutions */
#if SCIP_DISABLED_CODE
static
SCIP_DECL_CONSSEPASOL(consSepasolNonlinear)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of nonlinear constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consSepasolNonlinear NULL
#endif


/** constraint enforcing method of constraint handler for LP solutions */
static
SCIP_DECL_CONSENFOLP(consEnfolpNonlinear)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of nonlinear constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}


/** constraint enforcing method of constraint handler for relaxation solutions */
static
SCIP_DECL_CONSENFORELAX(consEnforelaxNonlinear)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of nonlinear constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}


/** constraint enforcing method of constraint handler for pseudo solutions */
static
SCIP_DECL_CONSENFOPS(consEnfopsNonlinear)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of nonlinear constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}


/** feasibility check method of constraint handler for integral solutions */
static
SCIP_DECL_CONSCHECK(consCheckNonlinear)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of nonlinear constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}


/** domain propagation method of constraint handler */
#if SCIP_DISABLED_CODE
static
SCIP_DECL_CONSPROP(consPropNonlinear)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of nonlinear constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consPropNonlinear NULL
#endif


/** presolving method of constraint handler */
#if SCIP_DISABLED_CODE
static
SCIP_DECL_CONSPRESOL(consPresolNonlinear)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of nonlinear constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consPresolNonlinear NULL
#endif


/** propagation conflict resolving method of constraint handler */
#if SCIP_DISABLED_CODE
static
SCIP_DECL_CONSRESPROP(consRespropNonlinear)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of nonlinear constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consRespropNonlinear NULL
#endif


/** variable rounding lock method of constraint handler */
static
SCIP_DECL_CONSLOCK(consLockNonlinear)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of nonlinear constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}


/** constraint activation notification method of constraint handler */
#if SCIP_DISABLED_CODE
static
SCIP_DECL_CONSACTIVE(consActiveNonlinear)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of nonlinear constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consActiveNonlinear NULL
#endif


/** constraint deactivation notification method of constraint handler */
#if SCIP_DISABLED_CODE
static
SCIP_DECL_CONSDEACTIVE(consDeactiveNonlinear)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of nonlinear constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consDeactiveNonlinear NULL
#endif


/** constraint enabling notification method of constraint handler */
#if SCIP_DISABLED_CODE
static
SCIP_DECL_CONSENABLE(consEnableNonlinear)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of nonlinear constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consEnableNonlinear NULL
#endif


/** constraint disabling notification method of constraint handler */
#if SCIP_DISABLED_CODE
static
SCIP_DECL_CONSDISABLE(consDisableNonlinear)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of nonlinear constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consDisableNonlinear NULL
#endif

/** variable deletion of constraint handler */
#if SCIP_DISABLED_CODE
static
SCIP_DECL_CONSDELVARS(consDelvarsNonlinear)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of nonlinear constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consDelvarsNonlinear NULL
#endif


/** constraint display method of constraint handler */
#if SCIP_DISABLED_CODE
static
SCIP_DECL_CONSPRINT(consPrintNonlinear)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of nonlinear constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consPrintNonlinear NULL
#endif


/** constraint copying method of constraint handler */
#if SCIP_DISABLED_CODE
static
SCIP_DECL_CONSCOPY(consCopyNonlinear)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of nonlinear constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consCopyNonlinear NULL
#endif


/** constraint parsing method of constraint handler */
#if SCIP_DISABLED_CODE
static
SCIP_DECL_CONSPARSE(consParseNonlinear)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of nonlinear constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consParseNonlinear NULL
#endif


/** constraint method of constraint handler which returns the variables (if possible) */
#if SCIP_DISABLED_CODE
static
SCIP_DECL_CONSGETVARS(consGetVarsNonlinear)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of nonlinear constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consGetVarsNonlinear NULL
#endif

/** constraint method of constraint handler which returns the number of variables (if possible) */
#if SCIP_DISABLED_CODE
static
SCIP_DECL_CONSGETNVARS(consGetNVarsNonlinear)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of nonlinear constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consGetNVarsNonlinear NULL
#endif

/** constraint handler method to suggest dive bound changes during the generic diving algorithm */
#if SCIP_DISABLED_CODE
static
SCIP_DECL_CONSGETDIVEBDCHGS(consGetDiveBdChgsNonlinear)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of nonlinear constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consGetDiveBdChgsNonlinear NULL
#endif

/** output method of statistics table to output file stream 'file' */
static
SCIP_DECL_TABLEOUTPUT(tableOutputNonlinear)
{ /*lint --e{715}*/
   SCIP_CONSHDLR* conshdlr;

   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);
   assert(conshdlr != NULL);

   /* print statistics for nonlinear handlers */
   printNlhdlrStatistics(scip, conshdlr, file);

   /* print statistics for constraint handler */
   printConshdlrStatistics(scip, conshdlr, file);

   return SCIP_OKAY;
}

/*
 * constraint specific interface methods
 */

/** creates the handler for nonlinear constraints and includes it in SCIP */
SCIP_RETCODE SCIPincludeConshdlrNonlinear(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata;

   /* create nonlinear constraint handler data */
   SCIP_CALL( SCIPallocClearMemory(scip, &conshdlrdata) );
   conshdlrdata->intevalvar = intEvalVarBoundTightening;
   conshdlrdata->curboundstag = 1;
   conshdlrdata->lastboundrelax = 1;
   conshdlrdata->curpropboundstag = 1;
   SCIP_CALL( SCIPcreateClock(scip, &conshdlrdata->canonicalizetime) );
   SCIP_CALL( SCIPqueueCreate(&conshdlrdata->reversepropqueue, 100, 2.0) );
   SCIP_CALL( SCIPhashmapCreate(&conshdlrdata->var2expr, SCIPblkmem(scip), 100) );

   /* include constraint handler */
   SCIP_CALL( SCIPincludeConshdlr(scip, CONSHDLR_NAME, CONSHDLR_DESC,
         CONSHDLR_SEPAPRIORITY, CONSHDLR_ENFOPRIORITY, CONSHDLR_CHECKPRIORITY,
         CONSHDLR_SEPAFREQ, CONSHDLR_PROPFREQ, CONSHDLR_EAGERFREQ, CONSHDLR_MAXPREROUNDS,
         CONSHDLR_DELAYSEPA, CONSHDLR_DELAYPROP, CONSHDLR_NEEDSCONS,
         CONSHDLR_PROP_TIMING, CONSHDLR_PRESOLTIMING,
         conshdlrCopyNonlinear,
         consFreeNonlinear, consInitNonlinear, consExitNonlinear,
         consInitpreNonlinear, consExitpreNonlinear, consInitsolNonlinear, consExitsolNonlinear,
         consDeleteNonlinear, consTransNonlinear, consInitlpNonlinear,
         consSepalpNonlinear, consSepasolNonlinear, consEnfolpNonlinear, consEnforelaxNonlinear, consEnfopsNonlinear, consCheckNonlinear,
         consPropNonlinear, consPresolNonlinear, consRespropNonlinear, consLockNonlinear,
         consActiveNonlinear, consDeactiveNonlinear,
         consEnableNonlinear, consDisableNonlinear, consDelvarsNonlinear,
         consPrintNonlinear, consCopyNonlinear, consParseNonlinear,
         consGetVarsNonlinear, consGetNVarsNonlinear, consGetDiveBdChgsNonlinear, conshdlrdata) );

   /* add nonlinear constraint handler parameters */
   /* TODO organize into more subcategories */
   SCIP_CALL( SCIPaddIntParam(scip, "constraints/" CONSHDLR_NAME "/maxproprounds",
         "limit on number of propagation rounds for a set of constraints within one round of SCIP propagation",
         &conshdlrdata->maxproprounds, FALSE, 10, 0, INT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "constraints/" CONSHDLR_NAME "/propauxvars",
         "whether to check bounds of all auxiliary variable to seed reverse propagation",
         &conshdlrdata->propauxvars, TRUE, TRUE, NULL, NULL) );

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

   SCIP_CALL( SCIPaddIntParam(scip, "constraints/" CONSHDLR_NAME "/bilinmaxnauxexprs",
           "maximal number of auxiliary expressions per bilinear term",
           &conshdlrdata->bilinmaxnauxexprs, FALSE, BILIN_MAXNAUXEXPRS, 0, INT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "constraints/" CONSHDLR_NAME "/reformbinprods",
         "whether to reformulate products of binary variables during presolving",
         &conshdlrdata->reformbinprods, FALSE, TRUE, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "constraints/" CONSHDLR_NAME "/reformbinprodsand",
         "whether to use the AND constraint handler for reformulating binary products",
         &conshdlrdata->reformbinprodsand, FALSE, TRUE, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(scip, "constraints/" CONSHDLR_NAME "/reformbinprodsfac",
         "minimum number of terms to reformulate bilinear binary products by factorizing variables (<= 1: disabled)",
         &conshdlrdata->reformbinprodsfac, FALSE, 50, 1, INT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "constraints/" CONSHDLR_NAME "/forbidmultaggrnlvar",
         "whether to forbid multiaggregation of nonlinear variables",
         &conshdlrdata->forbidmultaggrnlvar, TRUE, TRUE, NULL, NULL) );

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
         "signals a bound change to a nonlinear constraint", processVarEvent, NULL) );
   assert(conshdlrdata->eventhdlr != NULL);

   /* include table for statistics */
   assert(SCIPfindTable(scip, TABLE_NAME_NONLINEAR) == NULL);
   SCIP_CALL( SCIPincludeTable(scip, TABLE_NAME_NONLINEAR, TABLE_DESC_NONLINEAR, TRUE,
         NULL, NULL, NULL, NULL, NULL, NULL, tableOutputNonlinear,
         NULL, TABLE_POSITION_NONLINEAR, TABLE_EARLIEST_STAGE_NONLINEAR) );

   return SCIP_OKAY;
}

/** includes an expression constraint upgrade method into the expression constraint handler */
SCIP_RETCODE SCIPincludeConsUpgradeNonlinear(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_DECL_NONLINCONSUPGD((*nlconsupgd)),  /**< method to call for upgrading nonlinear constraint */
   int                   priority,           /**< priority of upgrading method */
   SCIP_Bool             active,             /**< should the upgrading method by active by default? */
   const char*           conshdlrname        /**< name of the constraint handler */
   )
{
   SCIP_CONSHDLR*     conshdlr;
   SCIP_CONSHDLRDATA* conshdlrdata;
   CONSUPGRADE*       consupgrade;
   char               paramname[SCIP_MAXSTRLEN];
   char               paramdesc[SCIP_MAXSTRLEN];
   int                i;

   assert(conshdlrname != NULL );
   assert(nlconsupgd != NULL);

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
   for( i = conshdlrdata->nconsupgrades - 1; i >= 0; --i )
   {
      if( conshdlrdata->consupgrades[i]->consupgd == nlconsupgd )
      {
#ifdef SCIP_DEBUG
         SCIPwarningMessage(scip, "Try to add already known upgrade method %p for constraint handler <%s>.\n", (void*)nlconsupgd, conshdlrname);
#endif
         return SCIP_OKAY;
      }
   }

   /* create an expression constraint upgrade data object */
   SCIP_CALL( SCIPallocBlockMemory(scip, &consupgrade) );
   consupgrade->consupgd = nlconsupgd;
   consupgrade->priority = priority;
   consupgrade->active   = active;

   /* insert expression constraint upgrade method into constraint handler data */
   SCIPensureBlockMemoryArray(scip, &conshdlrdata->consupgrades, &conshdlrdata->consupgradessize, conshdlrdata->nconsupgrades+1);
   assert(conshdlrdata->nconsupgrades+1 <= conshdlrdata->consupgradessize);

   for( i = conshdlrdata->nconsupgrades; i > 0 && conshdlrdata->consupgrades[i-1]->priority < consupgrade->priority; --i )
      conshdlrdata->consupgrades[i] = conshdlrdata->consupgrades[i-1];
   assert(0 <= i && i <= conshdlrdata->nconsupgrades);
   conshdlrdata->consupgrades[i] = consupgrade;
   conshdlrdata->nconsupgrades++;

   /* adds parameter to turn on and off the upgrading step */
   (void) SCIPsnprintf(paramname, SCIP_MAXSTRLEN, "constraints/" CONSHDLR_NAME "/upgrade/%s", conshdlrname);
   (void) SCIPsnprintf(paramdesc, SCIP_MAXSTRLEN, "enable nonlinear upgrading for constraint handler <%s>", conshdlrname);
   SCIP_CALL( SCIPaddBoolParam(scip,
         paramname, paramdesc,
         &consupgrade->active, FALSE, active, NULL, NULL) );

   return SCIP_OKAY;
}

/** creates and captures a nonlinear constraint
 *
 *  @note the constraint gets captured, hence at one point you have to release it using the method SCIPreleaseCons()
 */
SCIP_RETCODE SCIPcreateConsNonlinear(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS**           cons,               /**< pointer to hold the created constraint */
   const char*           name,               /**< name of constraint */
   SCIP_EXPR*            expr,               /**< expression of constraint (must not be NULL) */
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
   /* TODO: (optional) modify the definition of the SCIPcreateConsNonlinear() call, if you don't need all the information */
   SCIP_CONSHDLR* conshdlr;

   /* find the nonlinear constraint handler */
   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);
   if( conshdlr == NULL )
   {
      SCIPerrorMessage("nonlinear constraint handler not found\n");
      return SCIP_PLUGINNOTFOUND;
   }

   /* create constraint */
   SCIP_CALL( createCons(scip, conshdlr, cons, name, expr, lhs, rhs, TRUE,
      initial, separate, enforce, check, propagate, local, modifiable, dynamic, removable) );

   return SCIP_OKAY;
}

/** creates and captures a nonlinear constraint with all its constraint flags set to their
 *  default values
 *
 *  @note the constraint gets captured, hence at one point you have to release it using the method SCIPreleaseCons()
 */
SCIP_RETCODE SCIPcreateConsBasicNonlinear(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS**           cons,               /**< pointer to hold the created constraint */
   const char*           name,               /**< name of constraint */
   SCIP_EXPR*            expr,               /**< expression of constraint (must not be NULL) */
   SCIP_Real             lhs,                /**< left hand side of constraint */
   SCIP_Real             rhs                 /**< right hand side of constraint */
   )
{
   SCIP_CALL( SCIPcreateConsNonlinear(scip, cons, name, expr, lhs, rhs,
         TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   return SCIP_OKAY;
}

/** creates and captures a quadratic nonlinear constraint */
SCIP_RETCODE SCIPcreateConsQuadraticNonlinear(
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
   SCIP_EXPR* expr;

   assert(nlinvars == 0 || (linvars != NULL && lincoefs != NULL));
   assert(nquadterms == 0 || (quadvars1 != NULL && quadvars2 != NULL && quadcoefs != NULL));

   /* get expression constraint handler */
   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);
   if( conshdlr == NULL )
   {
      SCIPerrorMessage("nonlinear constraint handler not found\n");
      return SCIP_PLUGINNOTFOUND;
   }

   /* create quadratic expression */
   SCIP_CALL( SCIPcreateExprQuadratic(scip, &expr, nlinvars, linvars, lincoefs, nquadterms, quadvars1, quadvars2, quadcoefs, exprownerdataCreate, (SCIP_EXPR_OWNERDATACREATEDATA*)conshdlr) );
   assert(expr != NULL);

   /* create expression constraint */
   SCIP_CALL( createCons(scip, conshdlr, cons, name, expr, lhs, rhs, FALSE,
      initial, separate, enforce, check, propagate, local, modifiable, dynamic, removable) );

   /* release quadratic expression (captured by constraint now) */
   SCIP_CALL( SCIPreleaseExpr(scip, &expr) );

   return SCIP_OKAY;
}

/** returns the expression of the given nonlinear constraint */
SCIP_EXPR* SCIPgetExprConsNonlinear(
   SCIP_CONS*            cons                /**< constraint data */
   )
{
   SCIP_CONSDATA* consdata;

   assert(cons != NULL);
   assert(strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(cons)), CONSHDLR_NAME) == 0);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   return consdata->expr;
}

/** gets the left hand side of a nonlinear constraint */
SCIP_Real SCIPgetLhsConsNonlinear(
   SCIP_CONS*            cons                /**< constraint data */
   )
{
   SCIP_CONSDATA* consdata;

   assert(cons != NULL);
   assert(strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(cons)), CONSHDLR_NAME) == 0);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   return consdata->lhs;
}

/** gets the right hand side of a nonlinear constraint */
SCIP_Real SCIPgetRhsConsNonlinear(
   SCIP_CONS*            cons                /**< constraint data */
   )
{
   SCIP_CONSDATA* consdata;

   assert(cons != NULL);
   assert(strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(cons)), CONSHDLR_NAME) == 0);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   return consdata->rhs;
}

/** gets the nonlinear constraint as a nonlinear row representation. */
SCIP_RETCODE SCIPgetNlRowConsNonlinear(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint */
   SCIP_NLROW**          nlrow               /**< pointer to store nonlinear row */
   )
{
   SCIP_CONSDATA* consdata;

   assert(cons  != NULL);
   assert(nlrow != NULL);
   assert(strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(cons)), CONSHDLR_NAME) == 0);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   if( consdata->nlrow == NULL )
   {
//FIXME      SCIP_CALL( createNlRow(scip, cons) );
   }
   assert(consdata->nlrow != NULL);
   *nlrow = consdata->nlrow;

   return SCIP_OKAY;
}

/** returns the root curvature of the given nonlinear constraint
 *
 * @note The curvature information are computed during CONSINITSOL.
 */
SCIP_EXPRCURV SCIPgetCurvatureConsNonlinear(
   SCIP_CONS*            cons                /**< constraint data */
   )
{
   SCIP_CONSDATA* consdata;

   assert(cons != NULL);
   assert(strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(cons)), CONSHDLR_NAME) == 0);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   return consdata->curv;
}

/** returns representation of the expression of the given expression constraint as quadratic form, if possible
 *
 * Only sets *isquadratic to TRUE if the whole expression is quadratic (in the non-extended formulation) and non-linear.
 * That is, the expr in each SCIP_QUADEXPR_QUADTERM will be a variable expressions and
 * \ref SCIPgetVarExprVar() can be used to retrieve the variable.
 */
SCIP_RETCODE SCIPcheckQuadraticConsNonlinear(
   SCIP*                    scip,               /**< SCIP data structure */
   SCIP_CONS*               cons,               /**< constraint data */
   SCIP_Bool*               isquadratic         /**< buffer to store whether constraint is quadratic */
   )
{
   SCIP_CONSDATA* consdata;

   assert(scip != NULL);
   assert(cons != NULL);
   assert(isquadratic != NULL);
   assert(strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(cons)), CONSHDLR_NAME) == 0);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);
   assert(consdata->expr != NULL);

   /* check whether constraint expression is quadratic in extended formulation */
   SCIP_CALL( SCIPcheckExprQuadratic(scip, consdata->expr, isquadratic) );

   /* if not quadratic in non-extended formulation, then do indicate quadratic */
   if( *isquadratic )
      *isquadratic = SCIPexprAreQuadraticExprsVariables(consdata->expr);

   return SCIP_OKAY;
}

/** adds coef * var to expression constraint
 *
 * @attention This method can only be called in the problem stage.
 */
SCIP_RETCODE SCIPaddLinearTermConsNonlinear(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint data */
   SCIP_Real             coef,               /**< coefficient */
   SCIP_VAR*             var                 /**< variable */
   )
{
   SCIP_CONSHDLR* conshdlr;
   SCIP_CONSDATA* consdata;
   SCIP_EXPR* varexpr;

   assert(scip != NULL);
   assert(cons != NULL);

   if( SCIPgetStage(scip) != SCIP_STAGE_PROBLEM )
   {
      SCIPerrorMessage("SCIPaddLinearTermConsNonlinear can only be called in problem stage.\n");
      return SCIP_INVALIDCALL;
   }

   /* we should have an original constraint */
   assert(SCIPconsIsOriginal(cons));

   if( coef == 0.0 )
      return SCIP_OKAY;

   conshdlr = SCIPconsGetHdlr(cons);
   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);
   assert(consdata->expr != NULL);

   /* we should not have collected additional data for it
    * if some of these asserts fail, we may have to remove it and add some code to keep information uptodate
    */
   assert(consdata->nvarexprs == 0);
   assert(consdata->varexprs == NULL);
   assert(!consdata->catchedevents);

   SCIP_CALL( createExprVar(scip, conshdlr, &varexpr, var) );

   /* append to sum, if consdata->expr is sum and not used anywhere else */
   if( SCIPexprGetNUses(consdata->expr) == 1 && SCIPisExprSum(scip, consdata->expr) )
   {
      SCIP_CALL( SCIPappendExprSumExpr(scip, consdata->expr, varexpr, coef) );
   }
   else
   {
      /* create new expression = 1 * consdata->expr + coef * var */
      SCIP_EXPR* children[2] = { consdata->expr, varexpr };
      SCIP_Real coefs[2] = { 1.0, coef };

      SCIP_CALL( SCIPcreateExprSum(scip, &consdata->expr, 2, children, coefs, 0.0, exprownerdataCreate, (SCIP_EXPR_OWNERDATACREATEDATA*)conshdlr) );

      /* release old root expr */
      SCIP_CALL( SCIPreleaseExpr(scip, &children[0]) );
   }

   SCIP_CALL( SCIPreleaseExpr(scip, &varexpr) );

   /* not sure we care about any of these flags for original constraints */
   consdata->issimplified = FALSE;
   consdata->ispropagated = FALSE;

   return SCIP_OKAY;
}

/** returns an equivalent linear constraint if possible */
SCIP_RETCODE SCIPgetLinearConsNonlinear(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint data */
   SCIP_CONS**           lincons             /**< buffer to store linear constraint data */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_EXPR* expr;
   SCIP_VAR** vars;
   SCIP_Real lhs;
   SCIP_Real rhs;
   int i;

   assert(scip != NULL);
   assert(cons != NULL);
   assert(lincons != NULL);
   assert(strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(cons)), CONSHDLR_NAME) == 0);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);
   expr = consdata->expr;
   assert(expr != NULL);

   *lincons = NULL;

   /* not a linear constraint if the root expression is not a sum */
   if( !SCIPisExprSum(scip, expr) )
      return SCIP_OKAY;

   /* if at least one child is not a variable, then not a linear constraint */
   for( i = 0; i < SCIPexprGetNChildren(expr); ++i )
      if( !SCIPisExprVar(scip, SCIPexprGetChildren(expr)[i]) )
         return SCIP_OKAY;

   /* collect all variables */
   SCIP_CALL( SCIPallocBufferArray(scip, &vars, SCIPexprGetNChildren(expr)) );
   for( i = 0; i < SCIPexprGetNChildren(expr); ++i )
      vars[i] = SCIPgetVarExprVar(SCIPexprGetChildren(expr)[i]);

   /* consider constant part of the sum expression */
   lhs = SCIPisInfinity(scip, -consdata->lhs) ? -SCIPinfinity(scip) : (consdata->lhs - SCIPgetConstantExprSum(expr));
   rhs = SCIPisInfinity(scip,  consdata->rhs) ?  SCIPinfinity(scip) : (consdata->rhs - SCIPgetConstantExprSum(expr));

   SCIP_CALL( SCIPcreateConsLinear(scip, lincons, SCIPconsGetName(cons),
         SCIPexprGetNChildren(expr), vars, SCIPgetCoefsExprSum(expr),
         lhs, rhs,
         SCIPconsIsInitial(cons), SCIPconsIsSeparated(cons), SCIPconsIsEnforced(cons),
         SCIPconsIsChecked(cons), SCIPconsIsPropagated(cons), SCIPconsIsLocal(cons),
         SCIPconsIsModifiable(cons), SCIPconsIsDynamic(cons), SCIPconsIsRemovable(cons),
         SCIPconsIsStickingAtNode(cons)) );

   /* free memory */
   SCIPfreeBufferArray(scip, &vars);

   return SCIP_OKAY;
}
