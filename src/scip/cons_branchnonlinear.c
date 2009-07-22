//#define SCIP_DEBUG
#pragma ident "@(#) $Id: cons_branchnonlinear.c,v 1.1 2009/07/22 20:04:48 bzfviger Exp $"

/**@file    cons_branchnonlinear.c
 * @ingroup BRANCHINGRULES
 * @brief   constraint handler for branching on variables in nonlinear (nonconvex) constraints
 * @author  Stefan Vigerske
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>  /* for strcmp */

#include "scip/cons_linear.h"
#include "cons_branchnonlinear.h"


/* constraint handler properties */
#define CONSHDLR_NAME          "branchnonlinear"
#define CONSHDLR_DESC          "constraint handler for branching on nonlinear variables"
#define CONSHDLR_SEPAPRIORITY         0 /**< priority of the constraint handler for separation */
#define CONSHDLR_ENFOPRIORITY    -10000 /**< priority of the constraint handler for constraint enforcing: below cons_integral and nonlinear constraints */
#define CONSHDLR_CHECKPRIORITY        0 /**< priority of the constraint handler for checking feasibility */
#define CONSHDLR_SEPAFREQ            -1 /**< frequency for separating cuts; zero means to separate only in the root node */
#define CONSHDLR_PROPFREQ            -1 /**< frequency for propagating domains; zero means only preprocessing propagation */
#define CONSHDLR_EAGERFREQ           -1 /**< frequency for using all instead of only the useful constraints in separation,
                                              *   propagation and enforcement, -1 for no eager evaluations, 0 for first only */
#define CONSHDLR_MAXPREROUNDS         0 /**< maximal number of presolving rounds the constraint handler participates in (-1: no limit) */
#define CONSHDLR_DELAYSEPA        FALSE /**< should separation method be delayed, if other separators found cuts? */
#define CONSHDLR_DELAYPROP        FALSE /**< should propagation method be delayed, if other propagators found reductions? */
#define CONSHDLR_DELAYPRESOL      FALSE /**< should presolving method be delayed, if other presolvers found reductions? */
#define CONSHDLR_NEEDSCONS        FALSE /**< should the constraint handler be skipped, if no constraints are available? */

/*
 * Data structures
 */

struct VarInfeasibility
{
   SCIP_Real                 min;
   SCIP_Real                 max;
   SCIP_Real                 sum;
   struct VarInfeasibility*  next;
};
typedef struct VarInfeasibility VARINFEASIBILITY;

/** constraint handler data */
struct SCIP_ConshdlrData
{
   SCIP_HASHMAP*      branchcand;   /* branching candidates */
   VARINFEASIBILITY*  varinfeas;    /* list of variable infeasibilities */
   int                eventhdlrpos; /* filter position of event handler */
   
   char               strategy;     /* branching strategy */
   SCIP_Real          mindistbrpointtobound;  /* minimal (fractional) distance of branching point to bound */
};


/*
 * Local methods
 */

/* put your local methods here, and declare them static */

static
SCIP_RETCODE clearBranchingCandidates(SCIP* scip, SCIP_CONSHDLR* conshdlr)
{
   SCIP_CONSHDLRDATA* data;
   VARINFEASIBILITY* v;
   VARINFEASIBILITY* w = NULL;
  
   assert(scip != NULL);
   assert(conshdlr != NULL);
  
   data = SCIPconshdlrGetData(conshdlr);
   assert(data != NULL);
   assert(data->branchcand != NULL);
   
   v = data->varinfeas;
   while (v)
   {
      w = v->next;
      SCIPfreeBlockMemory(scip, &v);
      v = w;
   }
   data->varinfeas = NULL;
   
   SCIP_CALL( SCIPhashmapRemoveAll(data->branchcand) );
   
   return SCIP_OKAY;
}

static
SCIP_RETCODE selectBranchingPoint(SCIP* scip, SCIP_CONSHDLR* conshdlr, SCIP_VAR* var, SCIP_Real* leftub, SCIP_Real* rightlb)
{
   SCIP_CONSHDLRDATA*   data;
   SCIP_Real            branchpoint;
   SCIP_Real            lb, ub;

   assert(scip != NULL);
   assert(var  != NULL);
   assert(leftub  != NULL);
   assert(rightlb != NULL);
     
   assert(scip != NULL);
   assert(conshdlr != NULL);
  
   data = SCIPconshdlrGetData(conshdlr);
   assert(data != NULL);

   branchpoint = SCIPgetVarSol(scip, var);
   lb = SCIPvarGetLbLocal(var);
   ub = SCIPvarGetUbLocal(var);

   if (SCIPisInfinity(scip, branchpoint))
   {
      if (SCIPisPositive(scip, lb))
         branchpoint = lb + 1000;
      else
         branchpoint = 0.;
   }
   else if (SCIPisInfinity(scip, -branchpoint))
   {
      if (SCIPisNegative(scip, ub))
         branchpoint = ub - 1000;
      else
         branchpoint = 0.;
   }

   if (SCIPvarGetType(var) == SCIP_VARTYPE_CONTINUOUS)
   {
      if (!SCIPisInfinity(scip, -lb) && !SCIPisInfinity(scip, ub))
      {
         /* branch on value of LP solution
          * if it is too close to the bounds, move more into the middle of the interval */
         if (branchpoint < (1-data->mindistbrpointtobound) * lb + data->mindistbrpointtobound * ub)
            branchpoint = (1-data->mindistbrpointtobound) * lb + data->mindistbrpointtobound * ub;
         else if (branchpoint > data->mindistbrpointtobound * lb + (1-data->mindistbrpointtobound) * ub)
            branchpoint = data->mindistbrpointtobound * lb + (1-data->mindistbrpointtobound) * ub;

         /* for very tiny intervals we set it into the middle */
         if (!SCIPisGT(scip, branchpoint, lb) || !SCIPisLT(scip, branchpoint, ub))
            branchpoint = (lb+ub) * .5;
      }
      else if (!SCIPisLT(scip, lb, branchpoint))
      {
         assert(SCIPisInfinity(scip, ub));
         branchpoint = lb + MAX(0.5*ABS(lb), 1000);
      }
      else if (!SCIPisGT(scip, ub, branchpoint))
      {
         assert(SCIPisInfinity(scip, -lb));
         branchpoint = ub - MAX(0.5*ABS(ub), 1000);
      }

      *leftub = *rightlb = branchpoint;
   }
   else
   {
      if (branchpoint > ub)
         branchpoint = ub;
      else if (branchpoint < lb)
         branchpoint = lb;
      if (SCIPisIntegral(scip, branchpoint))
      {
         if (branchpoint < .5*(lb+ub))
            branchpoint += .5;
         else
            branchpoint -= .5;
      }
      *rightlb = SCIPceil(scip, branchpoint);
      *leftub  = SCIPfloor(scip, branchpoint);
   }
   
   return SCIP_OKAY;
}

static
SCIP_RETCODE selectBranchingVariable(SCIP* scip, SCIP_CONSHDLR* conshdlr, SCIP_VAR** var, SCIP_Real* leftub, SCIP_Real* rightlb)
{
   SCIP_CONSHDLRDATA*   data;
   int                  listidx;
   SCIP_HASHMAPLIST*    candlist;
   SCIP_VAR*            cand;
   VARINFEASIBILITY*    infeas;
   SCIP_Real            candleftub, candrightlb;
   SCIP_Real            deltaminus, deltaplus;
   SCIP_Real            pscostdown, pscostup;
   SCIP_Real            score, bestscore;
  
   assert(scip     != NULL);
   assert(conshdlr != NULL);
   assert(var      != NULL);
   assert(leftub   != NULL);
   assert(rightlb  != NULL);
  
   data = SCIPconshdlrGetData(conshdlr);
   assert(data != NULL);
   assert(data->branchcand != NULL);
   
   *var = NULL;
   bestscore = -1.0;
   
   for (listidx = 0; listidx < SCIPhashmapGetNLists(data->branchcand); ++listidx)
   {
      for (candlist = SCIPhashmapGetList(data->branchcand, listidx); candlist; candlist = SCIPhashmapListGetNext(candlist))
      {
         cand   = (SCIP_VAR*) SCIPhashmapListGetOrigin(candlist);
         infeas = (VARINFEASIBILITY*) SCIPhashmapListGetImage(candlist);

         switch (data->strategy)
         {
            case 'b':
               SCIP_CALL( selectBranchingPoint(scip, conshdlr, cand, &candleftub, &candrightlb) );
               assert(candleftub == candrightlb || SCIPvarGetType(cand) <= SCIP_VARTYPE_IMPLINT);
              
               if (SCIPisInfinity(scip, -SCIPvarGetLbLocal(cand)))
                  deltaminus = SCIPisInfinity(scip, infeas->max) ? SCIPinfinity(scip) : 0.1 * infeas->sum + 0.8 * infeas->min + 1.3 * infeas->max;
               else
                  deltaminus = candleftub - SCIPvarGetLbLocal(cand);
              
               if (SCIPisInfinity(scip,  SCIPvarGetUbLocal(cand)))
                  deltaplus  = SCIPisInfinity(scip, infeas->max) ? SCIPinfinity(scip) : 0.1 * infeas->sum + 0.8 * infeas->min + 1.3 * infeas->max;
               else
                  deltaplus  = SCIPvarGetUbLocal(cand) - candrightlb;
              
               break;
              
            case 'r':
               SCIP_CALL( selectBranchingPoint(scip, conshdlr, cand, &candleftub, &candrightlb) );
               assert(candleftub == candrightlb || SCIPvarGetType(cand) <= SCIP_VARTYPE_IMPLINT);
            
               if (SCIPisInfinity(scip, -SCIPvarGetLbLocal(cand)))
                  deltaplus  = SCIPisInfinity(scip, infeas->max) ? SCIPinfinity(scip) : 0.1 * infeas->sum + 0.8 * infeas->min + 1.3 * infeas->max;
               else
                  deltaplus  = candleftub - SCIPvarGetLbLocal(cand);
            
               if (SCIPisInfinity(scip,  SCIPvarGetUbLocal(cand)))
                  deltaminus = SCIPisInfinity(scip, infeas->max) ? SCIPinfinity(scip) : 0.1 * infeas->sum + 0.8 * infeas->min + 1.3 * infeas->max;
               else
                  deltaminus = SCIPvarGetUbLocal(cand) - candrightlb;
            
               break;

            case 'i':
               deltaminus = deltaplus = 0.1 * infeas->sum + 0.8 * infeas->min + 1.3 * infeas->max;
               break;

            default :
               SCIPerrorMessage("branching strategy %c unknown\n", data->strategy);
               return SCIP_ERROR;
         }

         if (SCIPisInfinity(scip, deltaminus) || SCIPisInfinity(scip, deltaplus))
            score = SCIPinfinity(scip);
         else
         {
            pscostdown = SCIPgetVarPseudocost(scip, cand, -deltaminus);
            pscostup   = SCIPgetVarPseudocost(scip, cand,  deltaplus);
            score      = SCIPgetBranchScore(scip, cand, pscostdown, pscostup);
         }
         SCIPdebugMessage("branching score variable %s = %g; \tinfeas = %g; \ttype=%d  bestscore=%g\n", SCIPvarGetName(cand), score, 0.1 * infeas->sum + 0.8 * infeas->min + 1.3 * infeas->max, SCIPvarGetType(cand), bestscore);

         if ( SCIPisSumGT(scip, score, bestscore) )
         {
            bestscore = score;
            *var      = cand;
            *leftub   = candleftub;
            *rightlb  = candrightlb;
         }
         else if ( SCIPisSumEQ(scip, score, bestscore) && !(SCIPisInfinity(scip, -SCIPvarGetLbLocal(*var)) && SCIPisInfinity(scip, SCIPvarGetUbLocal(*var))))
         { // if best candidate so far is bounded or unbounded at atmost one side, maybe take new candidate
            if ( (SCIPisInfinity(scip, -SCIPvarGetLbLocal(cand)) || SCIPisInfinity(scip, SCIPvarGetUbLocal(cand))) &&
                 (SCIPisInfinity(scip, -SCIPvarGetLbLocal(*var)) || SCIPisInfinity(scip, SCIPvarGetUbLocal(*var))) )
            { /* if both variables are unbounded but one of them is bounded on one side, take the one with the larger bound on this side (hope that this avoids branching on always the same variable) */
               if (SCIPvarGetUbLocal(cand) > SCIPvarGetUbLocal(*var) ||
                  SCIPvarGetLbLocal(cand) < SCIPvarGetLbLocal(*var))
               {
                  *var = cand;
                  *leftub   = candleftub;
                  *rightlb  = candrightlb;
               }
            }
            else if (SCIPvarGetType(*var) == SCIPvarGetType(cand))
            { /* if both have the same type, take the one with larger diameter */
               if (SCIPisLT(scip, SCIPvarGetUbLocal(*var) - SCIPvarGetLbLocal(*var), SCIPvarGetUbLocal(cand) - SCIPvarGetLbLocal(cand)))
               {
                  *var = cand;
                  *leftub   = candleftub;
                  *rightlb  = candrightlb;
               }
            }
            else if (SCIPvarGetType(*var) > SCIPvarGetType(cand))
            { /* take the one with better type ("more discrete") */
               *var = cand;
               *leftub   = candleftub;
               *rightlb  = candrightlb;
            }
         }
      }
   }

   return SCIP_OKAY;
}

static
SCIP_DECL_EVENTEXEC(processNodeFocusedEvent)
{
   assert(scip != NULL);
   assert(event != NULL);
   assert(eventdata != NULL);
   assert(SCIPeventGetType(event) & SCIP_EVENTTYPE_NODEFOCUSED);
   assert(strcmp(SCIPconshdlrGetName((SCIP_CONSHDLR*)eventdata), CONSHDLR_NAME) == 0);
   
   SCIP_CALL( clearBranchingCandidates(scip, (SCIP_CONSHDLR*)eventdata) );
   
   return SCIP_OKAY;
}

/*
 * Callback methods of constraint rule
 */


/** destructor of constraint handler to free constraint handler data (called when SCIP is exiting) */
#if 1
static
SCIP_DECL_CONSFREE(consFreeBranchNonlinear)
{
   SCIP_CONSHDLRDATA* data;
   
   assert(scip != NULL);
   assert(conshdlr != NULL);
   
   data = SCIPconshdlrGetData(conshdlr);
   if (!data)
      return SCIP_OKAY;
   
   assert(data->branchcand == NULL);
   assert(data->varinfeas  == NULL);
   
   SCIPfreeMemory(scip, &data);
   
   return SCIP_OKAY;   
}
#else
#define consFreeBranchNonlinear NULL
#endif


/** initialization method of constraint handler (called after problem was transformed) */
#if 0
static
SCIP_DECL_CONSINIT(consInitBranchNonlinear)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of branch nonlinear constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consInitBranchNonlinear NULL
#endif


/** deinitialization method of constraint handler (called before transformed problem is freed) */
#if 0
static
SCIP_DECL_CONSEXIT(consExitBranchNonlinear)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of branch nonlinear constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consExitBranchNonlinear NULL
#endif


/** presolving initialization method of constraint handler (called when presolving is about to begin) */
#if 0
static
SCIP_DECL_CONSINITPRE(consInitpreBranchNonlinear)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of branch nonlinear constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consInitpreBranchNonlinear NULL
#endif


/** presolving deinitialization method of constraint handler (called after presolving has been finished) */
#if 0
static
SCIP_DECL_CONSEXITPRE(consExitpreBranchNonlinear)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of branch nonlinear constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consExitpreBranchNonlinear NULL
#endif


/** solving process initialization method of constraint handler (called when branch and bound process is about to begin) */
#if 1
static
SCIP_DECL_CONSINITSOL(consInitsolBranchNonlinear)
{
   SCIP_CONSHDLRDATA* data;
   SCIP_EVENTHDLR* eventhdlr;
   
   assert(scip != NULL);
   assert(conshdlr != NULL);
   
   data = SCIPconshdlrGetData(conshdlr);
   assert(data != NULL);
   assert(data->branchcand == NULL);
   
   /* TODO: what is a good estimate for the hashmap size? should the constraint handler notify about the number of potential candidates? */
   SCIP_CALL( SCIPhashmapCreate(&data->branchcand, SCIPblkmem(scip), MAX(1, SCIPgetNVars(scip))) );
   
   eventhdlr = SCIPfindEventhdlr(scip, CONSHDLR_NAME);
   assert(eventhdlr != NULL);
   SCIP_CALL( SCIPcatchEvent(scip, SCIP_EVENTTYPE_NODEFOCUSED, eventhdlr, (SCIP_EVENTDATA*)conshdlr, &data->eventhdlrpos) );

   return SCIP_OKAY;
}
#else
#define consInitsolBranchNonlinear NULL
#endif

/** solving process deinitialization method of constraint handler (called before branch and bound process data is freed) */
#if 1
static
SCIP_DECL_CONSEXITSOL(consExitsolBranchNonlinear)
{
   SCIP_CONSHDLRDATA* data;
   SCIP_EVENTHDLR*      eventhdlr;
   
   assert(scip != NULL);
   assert(conshdlr != NULL);
   
   data = SCIPconshdlrGetData(conshdlr);
   assert(data != NULL);

   SCIP_CALL( clearBranchingCandidates(scip, conshdlr) );
   if (data->branchcand != NULL)
      SCIPhashmapFree(&data->branchcand);
   
   eventhdlr = SCIPfindEventhdlr(scip, CONSHDLR_NAME);
   assert(eventhdlr != NULL);
   SCIP_CALL( SCIPdropEvent(scip, SCIP_EVENTTYPE_NODEFOCUSED, eventhdlr, (SCIP_EVENTDATA*)conshdlr, data->eventhdlrpos) );
   
   return SCIP_OKAY;
}
#else
#define consExitsolBranchNonlinear NULL
#endif

/** frees specific constraint data */
#if 0
static
SCIP_DECL_CONSDELETE(consDeleteBranchNonlinear)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of branch nonlinear constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consDeleteBranchNonlinear NULL
#endif


/** transforms constraint data into data belonging to the transformed problem */ 
#if 0
static
SCIP_DECL_CONSTRANS(consTransBranchNonlinear)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of branch nonlinear constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consTransBranchNonlinear NULL
#endif


/** LP initialization method of constraint handler */
#if 0
static
SCIP_DECL_CONSINITLP(consInitlpBranchNonlinear)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of branch nonlinear constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consInitlpBranchNonlinear NULL
#endif


/** separation method of constraint handler for LP solutions */
#if 0
static
SCIP_DECL_CONSSEPALP(consSepalpBranchNonlinear)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of branch nonlinear constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consSepalpBranchNonlinear NULL
#endif


/** separation method of constraint handler for arbitrary primal solutions */
#if 0
static
SCIP_DECL_CONSSEPASOL(consSepasolBranchNonlinear)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of branch nonlinear constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consSepasolBranchNonlinear NULL
#endif


/** constraint enforcing method of constraint handler for LP solutions */
static
SCIP_DECL_CONSENFOLP(consEnfolpBranchNonlinear)
{
   SCIP_CONSHDLRDATA*   data;
   SCIP_VAR*            brvar = NULL;
   SCIP_Real            leftub=0., rightlb=0.;
   SCIP_Real            leftobjest, rightobjest; 
   
   assert(scip != NULL);
   assert(conshdlr != NULL);
   
   data = SCIPconshdlrGetData(conshdlr);
   assert(data != NULL);
   
   *result = SCIP_FEASIBLE;
   if (data->varinfeas == NULL)
   { /* have no candidates for branching */
      return SCIP_OKAY;
   }

   SCIP_CALL( selectBranchingVariable(scip, conshdlr, &brvar, &leftub, &rightlb) );
   SCIP_CALL( clearBranchingCandidates(scip, conshdlr) );
   
   if (!brvar)
   {
      printf("branching variable selection failed to select a variable\n");
      return SCIP_OKAY;
   }
   
   leftobjest = SCIPcalcChildEstimate(scip, brvar, leftub);
   rightobjest = (leftub != rightlb) ? SCIPcalcChildEstimate(scip, brvar, rightlb) : leftobjest;
   if (leftobjest > SCIPinfinity(scip))
      leftobjest = SCIPinfinity(scip)/5.;
   if (rightobjest > SCIPinfinity(scip))
      rightobjest = leftobjest;

   if (SCIPvarGetStatus(brvar) == SCIP_VARSTATUS_MULTAGGR)
   {
      SCIP_NODE* node;
      SCIP_CONS* cons;
      SCIP_Real  val = 1.0;
      SCIPdebugMessage("branching on multiaggregated variable %s: new intervals: [%g, %g] [%g, %g]\n", SCIPvarGetName(brvar), SCIPvarGetLbLocal(brvar), leftub, rightlb, SCIPvarGetUbLocal(brvar));

      SCIP_CALL( SCIPcreateChild(scip, &node, 0.0, leftobjest) );
      SCIP_CALL( SCIPcreateConsLinear(scip, &cons, "branch", 1, &brvar, &val, SCIPvarGetLbLocal(brvar), leftub, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE) );
      SCIP_CALL( SCIPaddConsNode(scip, node, cons, NULL) );
      SCIP_CALL( SCIPreleaseCons(scip, &cons) );

      SCIP_CALL( SCIPcreateChild(scip, &node, 0.0, rightobjest) );
      SCIP_CALL( SCIPcreateConsLinear(scip, &cons, "branch", 1, &brvar, &val, rightlb, SCIPvarGetUbLocal(brvar), TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE) );
      SCIP_CALL( SCIPaddConsNode(scip, node, cons, NULL) );
      SCIP_CALL( SCIPreleaseCons(scip, &cons) );
   }
   else
   {
      if (SCIPvarGetType(brvar) != SCIP_VARTYPE_CONTINUOUS)
      {
         SCIPdebugMessage("branching on discrete variable %s\n", SCIPvarGetName(brvar));
         SCIP_CALL( SCIPbranchVar(scip, brvar, NULL, NULL, NULL) );
      }
      else
      {
         SCIP_NODE* node;
         SCIPdebugMessage("branching on continuous variable %s: new intervals: [%g, %g] [%g, %g]\n", SCIPvarGetName(brvar), SCIPvarGetLbLocal(brvar), leftub, rightlb, SCIPvarGetUbLocal(brvar));

         SCIP_CALL( SCIPcreateChild(scip, &node, 0.0, leftobjest) );
         SCIP_CALL( SCIPchgVarUbNode(scip, node, brvar, leftub) );

         SCIP_CALL( SCIPcreateChild(scip, &node, 0.0, rightobjest) );
         SCIP_CALL( SCIPchgVarLbNode(scip, node, brvar, rightlb) );
      }
   }

   *result = SCIP_BRANCHED;

   return SCIP_OKAY;
}


/** constraint enforcing method of constraint handler for pseudo solutions */
static
SCIP_DECL_CONSENFOPS(consEnfopsBranchNonlinear)
{
   SCIP_CALL( consEnfolpBranchNonlinear(scip, conshdlr, conss, nconss, nusefulconss, solinfeasible, result) );
   
   return SCIP_OKAY;
}

/** feasibility check method of constraint handler for integral solutions */
static
SCIP_DECL_CONSCHECK(consCheckBranchNonlinear)
{
   assert(result != NULL);
   
   *result = SCIP_FEASIBLE;
   
   return SCIP_OKAY;
}


/** domain propagation method of constraint handler */
#if 0
static
SCIP_DECL_CONSPROP(consPropBranchNonlinear)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of branch nonlinear constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consPropBranchNonlinear NULL
#endif


/** presolving method of constraint handler */
#if 0
static
SCIP_DECL_CONSPRESOL(consPresolBranchNonlinear)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of branch nonlinear constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consPresolBranchNonlinear NULL
#endif


/** propagation conflict resolving method of constraint handler */
#if 0
static
SCIP_DECL_CONSRESPROP(consRespropBranchNonlinear)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of branch nonlinear constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consRespropBranchNonlinear NULL
#endif


/** variable rounding lock method of constraint handler */
static
SCIP_DECL_CONSLOCK(consLockBranchNonlinear)
{
   return SCIP_OKAY;
}


/** constraint activation notification method of constraint handler */
#if 0
static
SCIP_DECL_CONSACTIVE(consActiveBranchNonlinear)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of branch nonlinear constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consActiveBranchNonlinear NULL
#endif


/** constraint deactivation notification method of constraint handler */
#if 0
static
SCIP_DECL_CONSDEACTIVE(consDeactiveBranchNonlinear)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of branch nonlinear constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consDeactiveBranchNonlinear NULL
#endif


/** constraint enabling notification method of constraint handler */
#if 0
static
SCIP_DECL_CONSENABLE(consEnableBranchNonlinear)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of branch nonlinear constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consEnableBranchNonlinear NULL
#endif


/** constraint disabling notification method of constraint handler */
#if 0
static
SCIP_DECL_CONSDISABLE(consDisableBranchNonlinear)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of branch nonlinear constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consDisableBranchNonlinear NULL
#endif


/** constraint display method of constraint handler */
#if 0
static
SCIP_DECL_CONSPRINT(consPrintBranchNonlinear)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of branch nonlinear constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consPrintBranchNonlinear NULL
#endif


/** constraint copying method of constraint handler */
#if 0
static
SCIP_DECL_CONSCOPY(consCopyBranchNonlinear)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of branch nonlinear constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consCopyBranchNonlinear NULL
#endif


/** constraint parsing method of constraint handler */
#if 0
static
SCIP_DECL_CONSPARSE(consParseBranchNonlinear)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of branch nonlinear constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consParseBranchNonlinear NULL
#endif

/*
 * constraint specific interface methods
 */

/** creates the handler for branch nonlinear constraints and includes it in SCIP */
SCIP_RETCODE SCIPincludeConshdlrBranchNonlinear(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata;

   /* create branch nonlinear constraint handler data */
   SCIP_CALL( SCIPallocMemory(scip, &conshdlrdata) );
   conshdlrdata->branchcand   = NULL;
   conshdlrdata->varinfeas    = NULL;
   conshdlrdata->eventhdlrpos = -1;

   /* include constraint handler */
   SCIP_CALL( SCIPincludeConshdlr(scip, CONSHDLR_NAME, CONSHDLR_DESC,
         CONSHDLR_SEPAPRIORITY, CONSHDLR_ENFOPRIORITY, CONSHDLR_CHECKPRIORITY,
         CONSHDLR_SEPAFREQ, CONSHDLR_PROPFREQ, CONSHDLR_EAGERFREQ, CONSHDLR_MAXPREROUNDS, 
         CONSHDLR_DELAYSEPA, CONSHDLR_DELAYPROP, CONSHDLR_DELAYPRESOL, CONSHDLR_NEEDSCONS,
         consFreeBranchNonlinear, consInitBranchNonlinear, consExitBranchNonlinear, 
         consInitpreBranchNonlinear, consExitpreBranchNonlinear, consInitsolBranchNonlinear, consExitsolBranchNonlinear,
         consDeleteBranchNonlinear, consTransBranchNonlinear, consInitlpBranchNonlinear,
         consSepalpBranchNonlinear, consSepasolBranchNonlinear, consEnfolpBranchNonlinear, consEnfopsBranchNonlinear, consCheckBranchNonlinear, 
         consPropBranchNonlinear, consPresolBranchNonlinear, consRespropBranchNonlinear, consLockBranchNonlinear,
         consActiveBranchNonlinear, consDeactiveBranchNonlinear, 
         consEnableBranchNonlinear, consDisableBranchNonlinear,
         consPrintBranchNonlinear, consCopyBranchNonlinear, consParseBranchNonlinear,
         conshdlrdata) );


   /* add branch nonlinear constraint handler parameters */
   SCIP_CALL( SCIPaddCharParam(scip, "constraints/"CONSHDLR_NAME"/strategy", "strategy to use for selecting branching variable: b: rb-int-br, r: rb-int-br-rev, i: rb-inf", &conshdlrdata->strategy, FALSE, 'r', "bri", NULL, NULL) );
   SCIP_CALL( SCIPaddRealParam(scip, "constraints/"CONSHDLR_NAME"/mindistbrpointtobound", "minimal fractional distance of branching point to variable bounds; a value of 0.5 leads to branching always in the middle of a bounded domain", &conshdlrdata->mindistbrpointtobound, FALSE, 0.2, 0.0001, 0.5, NULL, NULL) );

   SCIP_CALL( SCIPincludeEventhdlr(scip, CONSHDLR_NAME, "handles node focused event for branching rule on nonlinear variables",
      NULL, NULL, NULL, NULL, NULL, NULL, processNodeFocusedEvent, NULL) );

   return SCIP_OKAY;
}


/** Updates or initializes the infeasibility of a variable.
 * If called the first time for some variable, then this variable is added to the list of branching candidates.
 */
SCIP_RETCODE SCIPconshdlrBranchNonlinearUpdateVarInfeasibility(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler */
   SCIP_VAR*             var,                /**< variable */
   SCIP_Real             varinfeasibility    /**< infeasibility of variable */
   )
{
   SCIP_CONSHDLRDATA* data;
   VARINFEASIBILITY*    varinfeas;
  
   SCIPdebugMessage("register infeasibility %g for variable %s  [%g, %g]\n", varinfeasibility, SCIPvarGetName(var), SCIPvarGetLbLocal(var), SCIPvarGetUbLocal(var));
   
   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);
   assert(varinfeasibility >= 0.0);
   assert(!SCIPisEQ(scip, SCIPvarGetLbLocal(var), SCIPvarGetUbLocal(var)));
  
   data = SCIPconshdlrGetData(conshdlr);
   assert(data != NULL);
   assert(data->branchcand != NULL);
   
   varinfeas = (VARINFEASIBILITY*)SCIPhashmapGetImage(data->branchcand, (void*)var);
   
   if (varinfeas == NULL)
   {
      SCIP_CALL( SCIPallocBlockMemory(scip, &varinfeas) );
      varinfeas->min  = varinfeasibility;
      varinfeas->max  = varinfeasibility;
      varinfeas->sum  = varinfeasibility;
      varinfeas->next = data->varinfeas;
      data->varinfeas = varinfeas;
      
      SCIP_CALL( SCIPhashmapInsert(data->branchcand, (void*)var, varinfeas) );
   }
   else
   {
      varinfeas->sum += varinfeasibility;
      if (varinfeasibility < varinfeas->min)
         varinfeas->min = varinfeasibility;
      if (varinfeasibility > varinfeas->max)
         varinfeas->max = varinfeasibility;
   }
   
   return SCIP_OKAY;
}
