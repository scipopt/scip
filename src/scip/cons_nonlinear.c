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
 * @brief  constraint handler for nonlinear constraints  specified by algebraic expressions
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
#include "scip/cons_and.h"
#include "scip/cons_bounddisjunction.h"
#include "scip/heur_subnlp.h"
#include "scip/heur_trysol.h"
#include "scip/debug.h"


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

/* put your local methods here, and declare them static */

/** callback that creates data that this conshdlr wants to store in an expression */
static
SCIP_DECL_EXPR_OWNERDATACREATE(exprownerdataCreate)
{
   assert(scip != NULL);
   assert(expr != NULL);
   assert(ownerdata != NULL);

   SCIP_CALL( SCIPallocClearBlockMemory(scip, ownerdata) );

   (*ownerdata)->nenfos = -1;

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

   /* expression should not be enforced anymore */
   assert((*ownerdata)->nenfos == 0);
   assert((*ownerdata)->auxvar == NULL);

   SCIPfreeBlockMemory(scip, ownerdata);

   return SCIP_OKAY;
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


/*
 * constraint specific interface methods
 */

/** creates the handler for nonlinear constraints and includes it in SCIP */
SCIP_RETCODE SCIPincludeConshdlrNonlinear(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONSHDLR* conshdlr;

   /* create nonlinear constraint handler data */
   conshdlrdata = NULL;
   /* TODO: (optional) create constraint handler specific data here */

   conshdlr = NULL;

   /* include constraint handler */
#if SCIP_DISABLED_CODE
   /* use SCIPincludeConshdlr() if you want to set all callbacks explicitly and realize (by getting compiler errors) when
    * new callbacks are added in future SCIP versions
    */
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
#else
   /* use SCIPincludeConshdlrBasic() plus setter functions if you want to set callbacks one-by-one and your code should
    * compile independent of new callbacks being added in future SCIP versions
    */
   SCIP_CALL( SCIPincludeConshdlrBasic(scip, &conshdlr, CONSHDLR_NAME, CONSHDLR_DESC,
         CONSHDLR_ENFOPRIORITY, CONSHDLR_CHECKPRIORITY, CONSHDLR_EAGERFREQ, CONSHDLR_NEEDSCONS,
         consEnfolpNonlinear, consEnfopsNonlinear, consCheckNonlinear, consLockNonlinear,
         conshdlrdata) );
   assert(conshdlr != NULL);

   /* set non-fundamental callbacks via specific setter functions */
   SCIP_CALL( SCIPsetConshdlrActive(scip, conshdlr, consActiveNonlinear) );
   SCIP_CALL( SCIPsetConshdlrCopy(scip, conshdlr, conshdlrCopyNonlinear, consCopyNonlinear) );
   SCIP_CALL( SCIPsetConshdlrDeactive(scip, conshdlr, consDeactiveNonlinear) );
   SCIP_CALL( SCIPsetConshdlrDelete(scip, conshdlr, consDeleteNonlinear) );
   SCIP_CALL( SCIPsetConshdlrDelvars(scip, conshdlr, consDelvarsNonlinear) );
   SCIP_CALL( SCIPsetConshdlrDisable(scip, conshdlr, consDisableNonlinear) );
   SCIP_CALL( SCIPsetConshdlrEnable(scip, conshdlr, consEnableNonlinear) );
   SCIP_CALL( SCIPsetConshdlrExit(scip, conshdlr, consExitNonlinear) );
   SCIP_CALL( SCIPsetConshdlrExitpre(scip, conshdlr, consExitpreNonlinear) );
   SCIP_CALL( SCIPsetConshdlrExitsol(scip, conshdlr, consExitsolNonlinear) );
   SCIP_CALL( SCIPsetConshdlrFree(scip, conshdlr, consFreeNonlinear) );
   SCIP_CALL( SCIPsetConshdlrGetDiveBdChgs(scip, conshdlr, consGetDiveBdChgsNonlinear) );
   SCIP_CALL( SCIPsetConshdlrGetVars(scip, conshdlr, consGetVarsNonlinear) );
   SCIP_CALL( SCIPsetConshdlrGetNVars(scip, conshdlr, consGetNVarsNonlinear) );
   SCIP_CALL( SCIPsetConshdlrInit(scip, conshdlr, consInitNonlinear) );
   SCIP_CALL( SCIPsetConshdlrInitpre(scip, conshdlr, consInitpreNonlinear) );
   SCIP_CALL( SCIPsetConshdlrInitsol(scip, conshdlr, consInitsolNonlinear) );
   SCIP_CALL( SCIPsetConshdlrInitlp(scip, conshdlr, consInitlpNonlinear) );
   SCIP_CALL( SCIPsetConshdlrParse(scip, conshdlr, consParseNonlinear) );
   SCIP_CALL( SCIPsetConshdlrPresol(scip, conshdlr, consPresolNonlinear, CONSHDLR_MAXPREROUNDS, CONSHDLR_PRESOLTIMING) );
   SCIP_CALL( SCIPsetConshdlrPrint(scip, conshdlr, consPrintNonlinear) );
   SCIP_CALL( SCIPsetConshdlrProp(scip, conshdlr, consPropNonlinear, CONSHDLR_PROPFREQ, CONSHDLR_DELAYPROP,
         CONSHDLR_PROP_TIMING) );
   SCIP_CALL( SCIPsetConshdlrResprop(scip, conshdlr, consRespropNonlinear) );
   SCIP_CALL( SCIPsetConshdlrSepa(scip, conshdlr, consSepalpNonlinear, consSepasolNonlinear, CONSHDLR_SEPAFREQ, CONSHDLR_SEPAPRIORITY, CONSHDLR_DELAYSEPA) );
   SCIP_CALL( SCIPsetConshdlrTrans(scip, conshdlr, consTransNonlinear) );
   SCIP_CALL( SCIPsetConshdlrEnforelax(scip, conshdlr, consEnforelaxNonlinear) );

#endif

#ifdef LINCONSUPGD_PRIORITY
   if( SCIPfindConshdlr(scip,"linear") != NULL )
   {
      /* include the linear constraint upgrade in the linear constraint handler */
      SCIP_CALL( SCIPincludeLinconsUpgrade(scip, linconsUpgdNonlinear, LINCONSUPGD_PRIORITY, CONSHDLR_NAME) );
   }
#endif

   /* add nonlinear constraint handler parameters */
   /* TODO: (optional) add constraint handler specific parameters with SCIPaddTypeParam() here */

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
   int                   nvars,              /**< number of variables in the constraint */
   SCIP_VAR**            vars,               /**< array with variables of constraint entries */
   SCIP_Real*            coefs,              /**< array with coefficients of constraint entries */
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
   SCIP_Bool             removable,          /**< should the relaxation be removed from the LP due to aging or cleanup?
                                              *   Usually set to FALSE. Set to TRUE for 'lazy constraints' and 'user cuts'. */
   SCIP_Bool             stickingatnode      /**< should the constraint always be kept at the node where it was added, even
                                              *   if it may be moved to a more global node?
                                              *   Usually set to FALSE. Set to TRUE to for constraints that represent node data. */
   )
{
   /* TODO: (optional) modify the definition of the SCIPcreateConsNonlinear() call, if you don't need all the information */

   SCIP_CONSHDLR* conshdlr;
   SCIP_CONSDATA* consdata;

   SCIPerrorMessage("method of nonlinear constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527} --e{715}*/

   /* find the nonlinear constraint handler */
   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);
   if( conshdlr == NULL )
   {
      SCIPerrorMessage("nonlinear constraint handler not found\n");
      return SCIP_PLUGINNOTFOUND;
   }

   /* create constraint data */
   consdata = NULL;
   /* TODO: create and store constraint specific data here */

   /* create constraint */
   SCIP_CALL( SCIPcreateCons(scip, cons, name, conshdlr, consdata, initial, separate, enforce, check, propagate,
         local, modifiable, dynamic, removable, stickingatnode) );

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
   int                   nvars,              /**< number of variables in the constraint */
   SCIP_VAR**            vars,               /**< array with variables of constraint entries */
   SCIP_Real*            coefs,              /**< array with coefficients of constraint entries */
   SCIP_Real             lhs,                /**< left hand side of constraint */
   SCIP_Real             rhs                 /**< right hand side of constraint */
   )
{
   SCIP_CALL( SCIPcreateConsNonlinear(scip, cons, name, nvars, vars, coefs, lhs, rhs,
         TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE) );

   return SCIP_OKAY;
}
