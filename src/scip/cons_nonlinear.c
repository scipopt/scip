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
#include "scip/cons_and.h"
#include "scip/cons_bounddisjunction.h"
#include "scip/heur_subnlp.h"
#include "scip/heur_trysol.h"
#include "nlpi/nlpi_ipopt.h"  /* for SCIPsolveLinearProb */
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
   SCIP_Longint          curboundstag;       /**< tag indicating current variable bounds */
   SCIP_Longint          lastboundrelax;     /**< tag when bounds where most recently relaxed */
   SCIP_Longint          lastvaractivitymethodchange; /**< tag when method used to evaluate activity of variables changed last */
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
//FIXME         SCIP_CALL( SCIPcallNlhdlrExitSepaNonlinear(scip, nlhdlr, expr, mydata->enfos[e]->nlhdlrexprdata) );
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
      SCIPexprSetActivity(expr, activity, conshdlrdata->curboundstag);
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

/** given three points, constructs coefficient of equation for hyperplane generated by these three points
 * Three points a, b, and c are given.
 * Computes coefficients alpha, beta, gamma, and delta, such that a, b, and c, satisfy
 * alpha * x1 + beta * x2 + gamma * x3 = delta and gamma >= 0.0.
 */
static
SCIP_RETCODE computeHyperplaneThreePoints(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             a1,                 /**< first coordinate of a */
   SCIP_Real             a2,                 /**< second coordinate of a */
   SCIP_Real             a3,                 /**< third coordinate of a */
   SCIP_Real             b1,                 /**< first coordinate of b */
   SCIP_Real             b2,                 /**< second coordinate of b */
   SCIP_Real             b3,                 /**< third coordinate of b */
   SCIP_Real             c1,                 /**< first coordinate of c */
   SCIP_Real             c2,                 /**< second coordinate of c */
   SCIP_Real             c3,                 /**< third coordinate of c */
   SCIP_Real*            alpha,              /**< coefficient of first coordinate */
   SCIP_Real*            beta,               /**< coefficient of second coordinate */
   SCIP_Real*            gamma_,             /**< coefficient of third coordinate */
   SCIP_Real*            delta               /**< constant right-hand side */
   )
{
   assert(scip != NULL);
   assert(alpha != NULL);
   assert(beta  != NULL);
   assert(gamma_ != NULL);
   assert(delta != NULL);

   *alpha  = -b3*c2 + a3*(-b2+c2) + a2*(b3-c3) + b2*c3;
   *beta   = -(-b3*c1 + a3*(-b1+c1) + a1*(b3-c3) + b1*c3);
   *gamma_ = -a2*b1 + a1*b2 + a2*c1 - b2*c1 - a1*c2 + b1*c2;
   *delta  = -a3*b2*c1 + a2*b3*c1 + a3*b1*c2 - a1*b3*c2 - a2*b1*c3 + a1*b2*c3;

   /* SCIPdebugMsg(scip, "alpha: %g beta: %g gamma: %g delta: %g\n", *alpha, *beta, *gamma_, *delta); */

   if( SCIPisInfinity(scip, REALABS(*gamma_ * a3)) ||
      SCIPisInfinity(scip, REALABS(*gamma_ * b3)) ||
      SCIPisInfinity(scip, REALABS(*gamma_ * c3)) )
   {
      SCIPdebugMsg(scip, "activity above SCIP infinity\n");
      *delta  = 0.0;
      *alpha  = 0.0;
      *beta   = 0.0;
      *gamma_ = 0.0;
      return SCIP_OKAY;
   }

   /* check if hyperplane contains all three points (necessary because of numerical troubles) */
   if( !SCIPisRelEQ(scip, *alpha * a1 + *beta * a2 - *delta, -*gamma_ * a3) ||
      !SCIPisRelEQ(scip, *alpha * b1 + *beta * b2 - *delta, -*gamma_ * b3) ||
      !SCIPisRelEQ(scip, *alpha * c1 + *beta * c2 - *delta, -*gamma_ * c3) )
   {
      SCIP_Real m[9];
      SCIP_Real rhs[3];
      SCIP_Real x[3];
      SCIP_Bool success;

      /*
      SCIPdebugMsg(scip, "a = (%g,%g,%g) hyperplane: %g rhs %g EQdelta: %d\n", a1, a2, a3, *alpha * a1 + *beta * a2 - *delta, -*gamma_ * a3, SCIPisRelEQ(scip, *alpha * a1 + *beta * a2 - *delta, -*gamma_ * a3));
      SCIPdebugMsg(scip, "b = (%g,%g,%g) hyperplane: %g rhs %g EQdelta: %d\n", b1, b2, b3, *alpha * b1 + *beta * b2 - *delta, -*gamma_ * b3, SCIPisRelEQ(scip, *alpha * b1 + *beta * b2 - *delta, -*gamma_ * b3));
      SCIPdebugMsg(scip, "c = (%g,%g,%g) hyperplane: %g rhs %g EQdelta: %d\n", c1, c2, c3, *alpha * c1 + *beta * c2 - *delta, -*gamma_ * c3, SCIPisRelEQ(scip, *alpha * c1 + *beta * c2 - *delta, -*gamma_ * c3));
      */

      /* initialize matrix column-wise */
      m[0] = a1;
      m[1] = b1;
      m[2] = c1;
      m[3] = a2;
      m[4] = b2;
      m[5] = c2;
      m[6] = a3;
      m[7] = b3;
      m[8] = c3;

      rhs[0] = 1.0;
      rhs[1] = 1.0;
      rhs[2] = 1.0;

      SCIPdebugMsg(scip, "numerical troubles - try to solve the linear system via an LU factorization\n");

      /* solve the linear problem */
      SCIP_CALL( SCIPsolveLinearProb(3, m, rhs, x, &success) );
      /* assert(success); */

      *delta  = rhs[0];
      *alpha  = x[0];
      *beta   = x[1];
      *gamma_ = x[2];

      /* set all coefficients to zero if one of the points is not contained in the hyperplane; this ensures that we do
       * not add a cut to SCIP and that all assertions are trivially fulfilled
       */
      if( !success || !SCIPisRelEQ(scip, *alpha * a1 + *beta * a2 - *delta, -*gamma_ * a3) ||
         !SCIPisRelEQ(scip, *alpha * b1 + *beta * b2 - *delta, -*gamma_ * b3) ||
         !SCIPisRelEQ(scip, *alpha * c1 + *beta * c2 - *delta, -*gamma_ * c3) ) /*lint !e774*/
      {
         SCIPdebugMsg(scip, "could not resolve numerical difficulties\n");
         *delta  = 0.0;
         *alpha  = 0.0;
         *beta   = 0.0;
         *gamma_ = 0.0;
      }
   }

   if( *gamma_ < 0.0 )
   {
      *alpha  = -*alpha;
      *beta   = -*beta;
      *gamma_ = -*gamma_;
      *delta  = -*delta;
   }

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
            SCIP_CALL( computeHyperplaneThreePoints(scip, p2[0], p2[1], p2val, p3[0], p3[1], p3val, p4[0], p4[1], p4val,
               &alpha, &beta, &gamma_, &delta) );
            /* if not underestimating in p1, then go to next candidate */
            if( alpha * p1[0] + beta * p1[1] + gamma_ * p1val - delta > 0.0 )
               continue;
            break;

         case 2 :
            /* get hyperplane through p1, p3, p4 */
            SCIP_CALL( computeHyperplaneThreePoints(scip, p1[0], p1[1], p1val, p3[0], p3[1], p3val, p4[0], p4[1], p4val,
               &alpha, &beta, &gamma_, &delta) );
            /* if not underestimating in p2, then go to next candidate */
            if( alpha * p2[0] + beta * p2[1] + gamma_ * p2val - delta > 0.0 )
               continue;
            break;

         case 3 :
            /* get hyperplane through p1, p2, p4 */
            SCIP_CALL( computeHyperplaneThreePoints(scip, p1[0], p1[1], p1val, p2[0], p2[1], p2val, p4[0], p4[1], p4val,
               &alpha, &beta, &gamma_, &delta) );
            /* if not underestimating in p3, then go to next candidate */
            if( alpha * p3[0] + beta * p3[1] + gamma_ * p3val - delta > 0.0 )
               continue;
            break;

         case 4 :
            /* get hyperplane through p1, p2, p3 */
            SCIP_CALL( computeHyperplaneThreePoints(scip, p1[0], p1[1], p1val, p2[0], p2[1], p2val, p3[0], p3[1], p3val,
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

/*
 * Callback methods of constraint handler
 */

/* TODO: Implement all necessary constraint handler methods. The methods with #if SCIP_DISABLED_CODE ... #else #define ... are optional */

/** copy method for constraint handler plugins (called when SCIP copies plugins) */
#if 1
static
SCIP_DECL_CONSHDLRCOPY(conshdlrCopyNonlinear)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of nonlinear constraint handler not implemented yet\n");
//   SCIPABORT(); /*lint --e{527}*/
   *valid = TRUE;

   return SCIP_OKAY;
}
#else
#define conshdlrCopyNonlinear NULL
#endif

/** destructor of constraint handler to free constraint handler data (called when SCIP is exiting) */
static
SCIP_DECL_CONSFREE(consFreeNonlinear)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;
   int i;

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   for( i = 0; i < conshdlrdata->nnlhdlrs; ++i )
   {
      SCIP_NLHDLR* nlhdlr;

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

      SCIPfreeMemory(scip, &nlhdlr->name);
      SCIPfreeMemoryNull(scip, &nlhdlr->desc);

      SCIPfreeMemory(scip, &nlhdlr);
   }

   SCIPfreeBlockMemoryArrayNull(scip, &conshdlrdata->nlhdlrs, conshdlrdata->nlhdlrssize);
   conshdlrdata->nlhdlrssize = 0;

   /* free upgrade functions */
   for( i = 0; i < conshdlrdata->nconsupgrades; ++i )
   {
      assert(conshdlrdata->consupgrades[i] != NULL);
      SCIPfreeBlockMemory(scip, &conshdlrdata->consupgrades[i]);
   }
   SCIPfreeBlockMemoryArrayNull(scip, &conshdlrdata->consupgrades, conshdlrdata->consupgradessize);

   SCIP_CALL( SCIPfreeClock(scip, &conshdlrdata->canonicalizetime) );

   SCIPqueueFree(&conshdlrdata->reversepropqueue);

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

   assert(conshdlrdata->branchrandnumgen == NULL);

   assert(SCIPhashmapGetNElements(conshdlrdata->var2expr) == 0);
   SCIPhashmapFree(&conshdlrdata->var2expr);

   SCIPfreeMemory(scip, &conshdlrdata);
   SCIPconshdlrSetData(conshdlr, NULL);

   return SCIP_OKAY;
}


/** initialization method of constraint handler (called after problem was transformed) */
#if 1
static
SCIP_DECL_CONSINIT(consInitNonlinear)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of nonlinear constraint handler not implemented yet\n");
//   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consInitNonlinear NULL
#endif


/** deinitialization method of constraint handler (called before transformed problem is freed) */
#if 1
static
SCIP_DECL_CONSEXIT(consExitNonlinear)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of nonlinear constraint handler not implemented yet\n");
//   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consExitNonlinear NULL
#endif


/** presolving initialization method of constraint handler (called when presolving is about to begin) */
#if 1
static
SCIP_DECL_CONSINITPRE(consInitpreNonlinear)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of nonlinear constraint handler not implemented yet\n");
//   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consInitpreNonlinear NULL
#endif


/** presolving deinitialization method of constraint handler (called after presolving has been finished) */
#if 1
static
SCIP_DECL_CONSEXITPRE(consExitpreNonlinear)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of nonlinear constraint handler not implemented yet\n");
//   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consExitpreNonlinear NULL
#endif


/** solving process initialization method of constraint handler (called when branch and bound process is about to begin) */
#if 1
static
SCIP_DECL_CONSINITSOL(consInitsolNonlinear)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of nonlinear constraint handler not implemented yet\n");
//   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consInitsolNonlinear NULL
#endif


/** solving process deinitialization method of constraint handler (called before branch and bound process data is freed) */
#if 1
static
SCIP_DECL_CONSEXITSOL(consExitsolNonlinear)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of nonlinear constraint handler not implemented yet\n");
//   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consExitsolNonlinear NULL
#endif


/** frees specific constraint data */
#if 1
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
#if 1
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
#if 1
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
#if 1
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
#if 1
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
#if 1
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
#if 1
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
#if 1
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
#if 1
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
#if 1
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
#if 1
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
#if 1
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

/** computes a facet of the convex or concave envelope of a vertex polyhedral function
 *
 * If \f$ f(x) \f$ is vertex-polyhedral, then \f$ g \f$ is a convex underestimator if and only if
 * \f$ g(v^i) \leq f(v^i), \forall i \f$, where \f$ \{ v^i \}_{i = 1}^{2^n} \subseteq \mathbb R^n \f$ are the vertices
 * of the domain of \f$ x \f$, \f$ [\ell,u] \f$. Hence, we can compute a linear underestimator by solving the following
 * LP (we don't necessarily get a facet of the convex envelope, see Technical detail below):
 *
 * \f{align*}{
 *              \max \, & \alpha^T x^* + \beta \\
 *     s.t. \; & \alpha^T v^i + \beta \le f(v^i), \, \forall i = 1, \ldots, 2^n
 * \f}
 *
 * In principle, one would need to update the LP whenever the domain changes. However, \f$ [\ell,u] = T([0, 1]^n) \f$,
 * where \f$ T \f$ is an affine linear invertible transformation given by \f$ T(y)_i = (u_i - \ell_i) y_i + \ell_i \f$.
 * Working with the change of variables \f$ x = T(y) \f$ allows us to keep the constraints of the LP, even if the domain
 * changes. Indeed, after the change of variables, the problem is: find an affine underestimator \f$ g \f$ such that \f$
 * g(T(y)) \le f(T(y)) \f$, for all \f$ y \in [0, 1]^n \f$. Now \f$ f(T(y)) \f$ is componentwise affine, but still
 * satisfies that \f$ g \f$ is a valid underestimator if and only if \f$ g(T(u)) \leq f(T(u)), \forall u \in \{0, 1\}^n
 * \f$. So we now look for \f$ \bar g(y) := g(T(y)) = g(((u_i - \ell_i) y_i + \ell_i)_i) = \bar \alpha^T y + \bar \beta
 * \f$, where \f$ \bar \alpha_i = (u_i - \ell_i) \alpha_i \f$ and \f$ \bar \beta = \sum_i \alpha_i \ell_i + \beta \f$. So
 * we find \f$ \bar g \f$ by solving the LP:
 *
 * \f{align*}{
 *              \max \, & \bar \alpha^T T^{-1}(x^*) + \bar \beta \\
 *     s.t. \; & \bar \alpha^T u + \bar \beta \le f(T(u)), \, \forall u \in \{0, 1\}^n
 * \f}
 *
 * and recover \f$ g \f$ by solving \f$ \bar \alpha_i = (u_i - \ell_i) \alpha_i, \bar \beta = \sum_i \alpha_i \ell_i +
 * \beta \f$. Notice that \f$ f(T(u^i)) = f(v^i) \f$ so the right hand side doesn't change after the change of variables.
 *
 * Furthermore, the LP has more constraints than variables, so we solve its dual:
 * \f{align*}{
 *              \min \, & \sum_i \lambda_i f(v^i) \\
 *     s.t. \; & \sum_i \lambda_i u^i = T^{-1}(x^*) \\
 *             & \sum_i \lambda_i = 1 \\
 *             & \forall i, \, \lambda_i \geq 0
 * \f}
 *
 * In case we look for an overestimate, we do exactly the same, but have to maximize in the dual LP instead
 * of minimize.
 *
 * #### Technical and implementation details
 * -# \f$ U \f$ has exponentially many variables, so we only apply this separator for \f$ n \leq 14 \f$.
 * -# If the bounds are not finite, there is no underestimator. Also, \f$ T^{-1}(x^*) \f$ must be in the domain,
 * otherwise the dual is infeasible.
 * -# After a facet is computed, we check whether it is a valid facet (i.e. we check \f$ \alpha^T v + \beta \le f(v) \f$
 *  for every vertex \f$ v \f$). If we find a violation of at most ADJUSTFACETFACTOR * SCIPlpfeastol, then we weaken \f$
 *  \beta \f$ by this amount, otherwise, we discard the cut.
 * -# If a variable is fixed within tolerances, we replace it with its value and compute the facet of the remaining
 * expression. Note that since we are checking the cut for validity, this will never produce wrong result.
 * -# If \f$ x^* \f$ is in the boundary of the domain, then the LP has infinitely many solutions, some of which might
 * have very bad numerical properties. For this reason, we perturb \f$ x^* \f$ to be in the interior of the region.
 * Furthermore, for some interior points, there might also be infinitely many solutions (e.g. for \f$ x y \f$ in \f$
 * [0,1]^2 \f$ any point \f$ (x^*, y^*) \f$ such that \f$ y^* = 1 - x^* \f$ has infinitely many solutions). For this
 * reason, we perturb any given \f$ x^* \f$. The idea is to try to get a facet of the convex/concave envelope. This only
 * happens when the solution has \f$ n + 1 \f$ non zero \f$ \lambda \f$'s (i.e. the primal has a unique solution).
 * -# We need to compute \f$ f(v^i) \f$ for every vertex of \f$ [\ell,u] \f$. A vertex is encoded by a number between 0
 * and \f$ 2^n - 1 \f$, via its binary representation (0 bit is lower bound, 1 bit is upper bound), so we can compute
 * all these values by iterating between 0 and \f$ 2^n - 1 \f$.
 * -# To check that the computed cut is valid we do the following: we use a gray code to loop over the vertices
 * of the box domain w.r.t. unfixed variables in order to evaluate the underestimator. To ensure the validity of the
 * underestimator, we check whether \f$ \alpha v^i + \beta \le f(v^i) \f$ for every vertex \f$ v^i \f$ and adjust
 * \f$ \beta \f$ if the maximal violation is small.
 *
 * @todo the solution is a facet if all variables of the primal have positive reduced costs (i.e. the solution is
 * unique). In the dual, this means that there are \f$ n + 1 \f$ variables with positive value. Can we use this or some
 * other information to handle any of both cases (point in the boundary or point in the intersection of polytopes
 * defining different pieces of the convex envelope)? In the case where the point is in the boundary, can we use that
 * information to maybe solve another to find a facet? How do the polytopes defining the pieces where the convex
 * envelope is linear looks like, i.e, given a point in the interior of a facet of the domain, does the midpoint of the
 * segment joining \f$ x^* \f$ with the center of the domain, always belongs to the interior of one of those polytopes?
 */
SCIP_RETCODE SCIPcomputeFacetVertexPolyhedralNonlinear(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< nonlinear constraint handler */
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
      SCIPwarningMessage(scip, "SCIPcomputeFacetVertexPolyhedralNonlinear() called with %d nonfixed variables. Must be between [1,%d].\n", nvars, SCIP_MAXVERTEXPOLYDIM);
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
