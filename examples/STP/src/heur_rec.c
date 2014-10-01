/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2014 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   heur_rec.c
 * @brief  rec primal heuristic
 * @author Daniel rehfeldt
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>
#include <stdio.h>
#include "scip/scip.h"
#include "scip/scipdefplugins.h"
#include "scip/cons_linear.h"
#include "heur_rec.h"
#include "scip/pub_misc.h"
#include "probdata_stp.h"

#define HEUR_NAME             "rec"
#define HEUR_DESC             "LNS heuristic fixing all variables corresponding to edges used in at least one of several selected solutions"
#define HEUR_DISPCHAR         'R'
#define HEUR_PRIORITY         100
#define HEUR_FREQ             1
#define HEUR_FREQOFS          0
#define HEUR_MAXDEPTH         -1
#define HEUR_TIMING           (SCIP_HEURTIMING_DURINGLPLOOP | SCIP_HEURTIMING_AFTERLPLOOP | SCIP_HEURTIMING_AFTERNODE)
#define HEUR_USESSUBSCIP      TRUE  /**< does the heuristic use a secondary SCIP instance? */

#define DEFAULT_MAXNODES      5000LL         /* maximum number of nodes to be regarded in the subproblem                   */
#define DEFAULT_MINIMPROVE    0.001          /* factor by which Rec should at least improve the incumbent       */
#define DEFAULT_MINNODES      50LL           /* minimum number of nodes to regard in the subproblem                   */
#define DEFAULT_MINFIXINGRATE 0.6            /* minimum percentage of variables to be fixed                           */
#define DEFAULT_NODESOFS      500LL          /* number of nodes added to the contingent of the total nodes            */
#define DEFAULT_NODESQUOT     0.1            /* subproblem nodes in relation to nodes of the original problem         */
#define DEFAULT_LPLIMFAC      2.0            /* factor by which the limit on the number of LP depends on the node limit */
#define DEFAULT_NUSEDSOLS     5              /* number of solutions that will be taken into account                   */
#define DEFAULT_NWAITINGNODES 200LL          /* number of nodes without incumbent change heuristic should wait        */
#define DEFAULT_NWAITINGSOLS  5              /* minimum number of new solutions before executing the heuristic again  */
#define DEFAULT_DONTWAITATROOT FALSE         /* should the nwaitingnodes parameter be ignored at the root node?       */
#define DEFAULT_USELPROWS     FALSE          /* should subproblem be created out of the rows in the LP rows,
                                              * otherwise, the copy constructors of the constraints handlers are used */
#define DEFAULT_COPYCUTS      FALSE           /* if DEFAULT_USELPROWS is FALSE, then should all active cuts from the
                                               * cutpool of the original scip be copied to constraints of the subscip
                                               */
#define DEFAULT_PERMUTE       FALSE          /* should the subproblem be permuted to increase diversification?        */
#define HASHSIZE_SOLS         11113          /* size of hash table for solution tuples in rec heuristic         */

/* event handler properties */
#define EVENTHDLR_NAME         "Rec"
#define EVENTHDLR_DESC         "LP event handler for rec heuristic"

/*
 * Data structures
 */

typedef struct SolTuple SOLTUPLE;

/** primal heuristic data */
struct SCIP_HeurData
{
   SCIP_SOL*             prevlastsol;        /**< worst solution taken into account during the previous run         */
   SCIP_SOL*             prevbestsol;        /**< best solution during the previous run                             */
   int                   prevnsols;          /**< number of all solutions during the previous run                   */

   SCIP_Longint          ncalls;
   SCIP_Longint          nlastsols;          /**< number of solutions available during the last run                 */
   SCIP_Longint          maxnodes;           /**< maximum number of nodes to regard in the subproblem               */
   SCIP_Longint          minnodes;           /**< minimum number of nodes to regard in the subproblem               */
   SCIP_Longint          nodesofs;           /**< number of nodes added to the contingent of the total nodes        */
   SCIP_Longint          usednodes;          /**< nodes already used by rec in earlier calls                  */
   SCIP_Real             nodesquot;          /**< subproblem nodes in relation to nodes of the original problem     */

   int                   nusedsols;          /**< number of solutions that will be taken into account               */
   int                   nselectedsols;      /**< number of solutions actually selected */
   int                   nwaitingsols;       /**< number of new solutions before executing the heuristic again      */
   SCIP_Longint          nwaitingnodes;      /**< number of nodes without incumbent change heuristic should wait    */
   unsigned int          nfailures;          /**< number of failures since last successful call                     */
   SCIP_Real             minfixingrate;      /**< minimum percentage of integer variables that have to be fixed     */
   SCIP_Real             minimprove;         /**< factor by which Rec should at least improve the incumbent   */
   SCIP_Real             nodelimit;          /**< the nodelimit employed in the current sub-SCIP, for the event handler*/
   SCIP_Real             lplimfac;           /**< factor by which the limit on the number of LP depends on the node limit */
   SCIP_Bool             dontwaitatroot;     /**< should the nwaitingnodes parameter be ignored at the root node?   */
   unsigned int          randseed;           /**< seed value for random number generator                            */
   SCIP_HASHTABLE*       hashtable;          /**< hashtable used to store the solution tuples already used          */
   SOLTUPLE*             lasttuple;          /**< last tuple of solutions created by rec                      */
   SCIP_Bool             uselprows;          /**< should subproblem be created out of the rows in the LP rows?      */
   SCIP_Bool             copycuts;           /**< if uselprows == FALSE, should all active cuts from cutpool be copied
                                              *   to constraints in subproblem?                                     */
   SCIP_Bool             permute;            /**< should the subproblem be permuted to increase diversification?    */
};

/** n-tuple of solutions and their hashkey */
struct SolTuple
{
   int*                  indices;            /**< sorted array of solution indices                                 */
   int                   size;               /**< size of the array (should be heurdata->nusedsols)                */
   unsigned int          key;                /**< hashkey of the tuple                                             */
   SOLTUPLE*             prev;               /**< previous solution tuple created                                  */
};


/*
 * Local methods
 */


/** gets the hash key of a solution tuple */
static
SCIP_DECL_HASHGETKEY(hashGetKeySols)
{  /*lint --e{715}*/
   return elem;
}


/** returns TRUE iff both solution tuples are identical */
static
SCIP_DECL_HASHKEYEQ(hashKeyEqSols)
{  /*lint --e{715}*/
   int i;
   int size;

   int* indices1;
   int* indices2;

   indices1 = ((SOLTUPLE*)key1)->indices;
   indices2 = ((SOLTUPLE*)key2)->indices;

   /* there should be two nonempty arrays of the same size */
   assert(indices1 != NULL);
   assert(indices2 != NULL);
   assert(((SOLTUPLE*)key1)->size == ((SOLTUPLE*)key2)->size);

   size  = ((SOLTUPLE*)key1)->size;

   /* compare arrays by components, return TRUE, iff equal */
   for( i = 0; i < size; i++ )
   {
      if( indices1[i] != indices2[i] )
         return FALSE;
   }

   return TRUE;
}


/** returns hashkey of a solution tuple */
static
SCIP_DECL_HASHKEYVAL(hashKeyValSols)
{  /*lint --e{715}*/
   return ((SOLTUPLE*)key)->key;
}


/** calculates a hash key for a given tuple of solution indices */
static
unsigned int calculateHashKey(
   int*                  indices,            /**< indices of solutions */
   int                   size                /**< number of solutions */
   )
{
   int i;
   unsigned int hashkey;

   /* hashkey should be (x1+1) * (x2+1) * ... * (xn+1) + x1 + x2 + ... + xn */
   hashkey = 1;
   for( i = 0; i < size; i++ )
      hashkey *= indices[i] + 1;
   for( i = 0; i < size; i++ )
      hashkey += indices[i];

   return hashkey;
}


/** insertion sort for a small int array */
static void sortArray(
   int*                  a,                  /**< array to be sorted */
   int                   size                /**< size of array */
   )
{
   int i;
   int j;
   int tmp;

   /* simple insertion sort algorithm */
   for( i = 1; i < size; i++ )
   {
      tmp = a[i];
      j = i-1;
      while( j >= 0 && a[j] > tmp )
      {
         a[j+1] = a[j]; /*lint !e679*/
         j = j-1;
      }
      a[j+1] = tmp; /*lint !e679*/
   }
}


/** creates a new tuple of solutions */
static
SCIP_RETCODE createSolTuple(
   SCIP*                 scip,               /**< original SCIP data structure */
   SOLTUPLE**            elem,               /**< tuple of solutions which should be created */
   int*                  indices,            /**< indices of solutions */
   int                   size,               /**< number of solutions */
   SCIP_HEURDATA*        heurdata            /**< primal heuristic data */
   )
{
   /* memory allocation */
   SCIP_CALL( SCIPallocBlockMemory(scip, elem) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(*elem)->indices, size) );
   BMScopyMemoryArray((*elem)->indices, indices, size);

   /* data input */
   sortArray(indices, size);
   (*elem)->size = size;
   (*elem)->key = calculateHashKey((*elem)->indices, (*elem)->size);
   (*elem)->prev = heurdata->lasttuple;

   /* update heurdata */
   heurdata->lasttuple = *elem;
   return SCIP_OKAY;
}


/** creates all variables of the subproblem */
static SCIP_RETCODE fixVars(
   SCIP*                 scip,               /**< original SCIP data structure */
   SCIP*                 subscip,            /**< SCIP data structure for the subproblemd */
   SCIP_VAR**            subvars,            /**< the variables of the subproblem */
   int*                  selection,          /**< pool of solutions rec will use */
   SCIP_HEURDATA*        heurdata,           /**< primal heuristic data */
   SCIP_Bool*            success             /**< pointer to store whether the problem was created successfully */
   )
{
   SCIP_VAR** vars;                          /* original scip variables                */
   SCIP_SOL** sols;                          /* pool of solutions                      */
   SCIP_Real fixingrate;                     /* percentage of variables that are fixed */
   SCIP_Bool fixable;
   SCIP_Real varsolval;
   SCIP_Real varrevsolval;

   int nvars;
   int i;
   int j;
   int fixingcounter;

   sols = SCIPgetSols(scip);
   assert(sols != NULL);

   /* get required data of the original problem */
   vars = SCIPprobdataGetEdgeVars(scip);
   nvars = SCIPprobdataGetNVars(scip);
   assert(nvars == SCIPprobdataGetNEdges(scip));

   fixingcounter = 0;

   /* create the variables of the subproblem */
   for( i = 0; i < nvars; i += 2 )
   {
      fixable = TRUE;

      /* check, whether the value of the variable is identical for each selected solution */
      for( j = 0; j < heurdata->nselectedsols; j++ )
      {
         varsolval = SCIPgetSolVal(scip, sols[selection[j]], vars[i]);
	 varrevsolval = SCIPgetSolVal(scip, sols[selection[j]], vars[i+1]);
	 if( SCIPisEQ(scip, varsolval, 1.0) || SCIPisEQ(scip, varrevsolval, 1.0) )
	 {
	    fixable = FALSE;
	    break;
	 }
      }

      /* original solval can be outside transformed global bounds */
      //fixable = fixable && SCIPvarGetLbGlobal(vars[i]) <= solval && solval <= SCIPvarGetUbGlobal(vars[i]);

      /* if the current variable is not contained in any selected solution, fix it to zero */
      if( fixable )
      {
         SCIP_CALL( SCIPchgVarLbGlobal(subscip, subvars[i], 0.0) );
         SCIP_CALL( SCIPchgVarUbGlobal(subscip, subvars[i], 0.0) );
	 SCIP_CALL( SCIPchgVarLbGlobal(subscip, subvars[i+1], 0.0) );
         SCIP_CALL( SCIPchgVarUbGlobal(subscip, subvars[i+1], 0.0) );
         fixingcounter += 2;
      }

   }

   fixingrate = (SCIP_Real)fixingcounter / (SCIP_Real)(MAX(nvars, 1));
   //printf("fixingrate %f \n\n", fixingrate);

   /* if all variables have been fixed or the amount of fixed variables is insufficient, abort immediately */
   *success = fixingcounter < nvars; //&& fixingrate >= heurdata->minfixingrate;
   /*
     if( *success )
     printf("fixed %d vars \n\n", fixingcounter);
     else
     printf("no success fixingcounter: %d, nvars: %d\n\n", fixingcounter, nvars);
   */
   return SCIP_OKAY;
}

/** creates the rows of the subproblem */
static
SCIP_RETCODE createRows(
   SCIP*                 scip,               /**< original SCIP data structure */
   SCIP*                 subscip,            /**< SCIP data structure for the subproblem */
   SCIP_VAR**            subvars             /**< the variables of the subproblem */
   )
{
   SCIP_ROW** rows;                          /* original scip rows                       */
   SCIP_CONS* cons;                          /* new constraint                           */
   SCIP_VAR** consvars;                      /* new constraint's variables               */
   SCIP_COL** cols;                          /* original row's columns                   */

   SCIP_Real constant;                       /* constant added to the row                */
   SCIP_Real lhs;                            /* left hand side of the row                */
   SCIP_Real rhs;                            /* left right side of the row               */
   SCIP_Real* vals;                          /* variables' coefficient values of the row */

   int nrows;
   int nnonz;
   int i;
   int j;

   /* get the rows and their number */
   SCIP_CALL( SCIPgetLPRowsData(scip, &rows, &nrows) );

   /* copy all rows to linear constraints */
   for( i = 0; i < nrows; i++ )
   {
      /* ignore rows that are only locally valid */
      if( SCIProwIsLocal(rows[i]) )
         continue;

      /* get the row's data */
      constant = SCIProwGetConstant(rows[i]);
      lhs = SCIProwGetLhs(rows[i]) - constant;
      rhs = SCIProwGetRhs(rows[i]) - constant;
      vals = SCIProwGetVals(rows[i]);
      nnonz = SCIProwGetNNonz(rows[i]);
      cols = SCIProwGetCols(rows[i]);

      assert(lhs <= rhs);

      /* allocate memory array to be filled with the corresponding subproblem variables */
      SCIP_CALL( SCIPallocBufferArray(scip, &consvars, nnonz) );
      for( j = 0; j < nnonz; j++ )
         consvars[j] = subvars[SCIPvarGetProbindex(SCIPcolGetVar(cols[j]))];

      /* create a new linear constraint and add it to the subproblem */
      SCIP_CALL( SCIPcreateConsLinear(subscip, &cons, SCIProwGetName(rows[i]), nnonz, consvars, vals, lhs, rhs,
            TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, TRUE, TRUE, FALSE) );
      SCIP_CALL( SCIPaddCons(subscip, cons) );
      SCIP_CALL( SCIPreleaseCons(subscip, &cons) );

      /* free temporary memory */
      SCIPfreeBufferArray(scip, &consvars);
   }

   return SCIP_OKAY;
}

/** creates a subproblem for subscip by fixing certain variables */
static
SCIP_RETCODE setupSubproblem(
   SCIP*                 scip,               /**< original SCIP data structure */
   SCIP*                 subscip,            /**< SCIP data structure for the subproblem */
   SCIP_VAR**            subvars,            /**< the variables of the subproblem */
   int*                  selection,          /**< pool of solutions rec will use */
   SCIP_HEURDATA*        heurdata,           /**< primal heuristic data */
   SCIP_Bool*            success             /**< pointer to store whether the problem was created successfully */
   )
{
   SCIP_SOL** sols;                          /* array of all solutions found so far         */
   SOLTUPLE* elem;
   SCIP_Real* soltimes;
   int i;
   int j;
   int min;
   int nselectedsols = 0;
   int nsols;                                /* number of all solutions found so far        */
   int ncalls;                               /* number of heuristic calls */
   int nusedsols;                            /* number of solutions to use in rec     */
   int* solselected;
   char consname[SCIP_MAXSTRLEN];

   /* get solution data */
   nsols = SCIPgetNSols(scip);
   sols = SCIPgetSols(scip);
   ncalls = heurdata->ncalls;
   nusedsols = heurdata->nusedsols;
   assert(nusedsols > 1);
   assert(nsols >= nusedsols);

   SCIP_CALL( SCIPallocBufferArray(scip, &solselected, nsols) );
   SCIP_CALL( SCIPallocBufferArray(scip, &soltimes, nsols) );
   for( i = 0; i < nsols; i++ )
      solselected[i] = FALSE;

   selection[nselectedsols++] = 0;
   solselected[0] = TRUE;
   //printf("solution: %d, found by: %s \n", SCIPsolGetIndex(sols[0]), SCIPheurGetName(SCIPsolGetHeur(sols[0])));
   for( i = 1; i < nsols && nselectedsols < nusedsols - 2; i++ )
   {
      if( strcmp(SCIPheurGetName(SCIPsolGetHeur(sols[i - 1])), "local") != 0 )
      {
         selection[nselectedsols++] = i;
	 solselected[i] = TRUE;
         // printf("solution: %d, found by: %s \n", SCIPsolGetIndex(sols[i]), SCIPheurGetName(SCIPsolGetHeur(sols[i])));
      }
   }

   for( i = 0; i < nsols; i++ )
      soltimes[i] = SCIPsolGetTime(sols[i]);

   SCIPsortRealIntPtr(soltimes, solselected, (void**) sols,  nsols);

   SCIPfreeBufferArray(scip, &soltimes);

   min = MIN(nusedsols, heurdata->nwaitingsols);

   for( i = 0; i < min && nselectedsols < nusedsols; i++ )
   {
      j = nsols - ((i + ncalls) % min) - 1;
      if( solselected[j] == FALSE )
      {
         selection[nselectedsols++] = j;
         solselected[j] = TRUE;
         //   printf("solution(phase min): %d found by: %s \n", SCIPsolGetIndex(sols[j]), SCIPheurGetName(SCIPsolGetHeur(sols[j])));
      }
   }

   for( i = 0; i < nsols && nselectedsols < nusedsols; i++ )
   {
      j = (i + ncalls + nusedsols) % nsols;
      if( solselected[j] == FALSE && strcmp(SCIPheurGetName(SCIPsolGetHeur(sols[j - 1])), "local") != 0 )
      {
         selection[nselectedsols++] = j;
         solselected[j] = TRUE;
         //    printf("solution(phase 2): %d found by: %s \n", j, SCIPheurGetName(SCIPsolGetHeur(sols[j])));
      }
   }
   heurdata->nselectedsols = nselectedsols;


   SCIPfreeBufferArray(scip, &solselected);
   SCIP_CALL( createSolTuple(scip, &elem, selection, nselectedsols, heurdata) );

   *success = !SCIPhashtableExists(heurdata->hashtable, elem);

   /* check, whether solution tuple has already been tried */
   if( !SCIPhashtableExists(heurdata->hashtable, elem) )
   {
      SCIP_CALL( SCIPhashtableInsert(heurdata->hashtable, elem) );
   }

   /* no acceptable solution tuple could be created */
   if( !(*success) )
      return SCIP_OKAY;

   /* get name of the original problem and add the string "_recsub" */
   (void) SCIPsnprintf(consname, SCIP_MAXSTRLEN, "%s_recsub", SCIPgetProbName(scip));

   /* set up the variables of the subproblem */
   SCIP_CALL( fixVars(scip, subscip, subvars, selection, heurdata, success) );

   /* we copy the rows of the LP, if the enough variables could be fixed and we work on the MIP
      relaxation of the problem */
   if( *success && heurdata->uselprows )
   {
      SCIP_CALL( createRows(scip, subscip, subvars) );
   }

   return SCIP_OKAY;
}


/** creates a new solution for the original problem by copying the solution of the subproblem */
static
SCIP_RETCODE createNewSol(
   SCIP*                 scip,               /**< original SCIP data structure */
   SCIP*                 subscip,            /**< SCIP structure of the subproblem */
   SCIP_VAR**            subvars,            /**< the variables of the subproblem */
   SCIP_HEUR*            heur,               /**< rec heuristic structure */
   SCIP_SOL*             subsol,             /**< solution of the subproblem */
   int*                  solindex,           /**< index of the solution */
   SCIP_Bool*            success             /**< used to store whether new solution was found or not */
   )
{
   SCIP_VAR** vars;                          /* the original problem's variables                */
   int        nvars;
   SCIP_SOL*  newsol;                        /* solution to be created for the original problem */
   SCIP_Real* subsolvals;                    /* solution values of the subproblem               */

   assert(scip != NULL);
   assert(subscip != NULL);
   assert(subvars != NULL);
   assert(subsol != NULL);

   /* get variables' data */
   SCIP_CALL( SCIPgetVarsData(scip, &vars, &nvars, NULL, NULL, NULL, NULL) );
   /* sub-SCIP may have more variables than the number of active (transformed) variables in the main SCIP
    * since constraint copying may have required the copy of variables that are fixed in the main SCIP
    */
   assert(nvars <= SCIPgetNOrigVars(subscip));

   SCIP_CALL( SCIPallocBufferArray(scip, &subsolvals, nvars) );

   /* copy the solution */
   SCIP_CALL( SCIPgetSolVals(subscip, subsol, nvars, subvars, subsolvals) );

   /* create new solution for the original problem */
   SCIP_CALL( SCIPcreateSol(scip, &newsol, heur) );
   SCIP_CALL( SCIPsetSolVals(scip, newsol, nvars, vars, subsolvals) );
   *solindex = SCIPsolGetIndex(newsol);

   /* try to add new solution to scip and free it immediately */
   SCIP_CALL( SCIPtrySolFree(scip, &newsol, FALSE, TRUE, TRUE, TRUE, success) );

   SCIPfreeBufferArray(scip, &subsolvals);

   return SCIP_OKAY;
}


/* ---------------- Callback methods of event handler ---------------- */

/* exec the event handler
 *
 * we interrupt the solution process
 */
static
SCIP_DECL_EVENTEXEC(eventExecRec)
{
   SCIP_HEURDATA* heurdata;

   assert(eventhdlr != NULL);
   assert(eventdata != NULL);
   assert(strcmp(SCIPeventhdlrGetName(eventhdlr), EVENTHDLR_NAME) == 0);
   assert(event != NULL);
   assert(SCIPeventGetType(event) & SCIP_EVENTTYPE_LPSOLVED);

   heurdata = (SCIP_HEURDATA*)eventdata;
   assert(heurdata != NULL);

   /* interrupt solution process of sub-SCIP */
   if( SCIPgetNLPs(scip) > heurdata->lplimfac * heurdata->nodelimit )
   {
      SCIPdebugMessage("interrupt after  %"SCIP_LONGINT_FORMAT" LPs\n",SCIPgetNLPs(scip));
      SCIP_CALL( SCIPinterruptSolve(scip) );
   }

   return SCIP_OKAY;
}

/*
 * Callback methods of primal heuristic
 */

/** copy method for primal heuristic plugins (called when SCIP copies plugins) */
static
SCIP_DECL_HEURCOPY(heurCopyRec)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(heur != NULL);
   assert(strcmp(SCIPheurGetName(heur), HEUR_NAME) == 0);

   /* call inclusion method of primal heuristic */
   SCIP_CALL( SCIPincludeHeurRec(scip) );

   return SCIP_OKAY;
}

/** destructor of primal heuristic to free user data (called when SCIP is exiting) */
static
SCIP_DECL_HEURFREE(heurFreeRec)
{  /*lint --e{715}*/
   SCIP_HEURDATA* heurdata;

   assert(heur != NULL);
   assert(scip != NULL);

   /* get heuristic data */
   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);

   /* free heuristic data */
   SCIPfreeMemory(scip, &heurdata);
   SCIPheurSetData(heur, NULL);

   return SCIP_OKAY;
}

/** initialization method of primal heuristic (called after problem was transformed) */
static
SCIP_DECL_HEURINIT(heurInitRec)
{  /*lint --e{715}*/
   SCIP_HEURDATA* heurdata;

   assert(heur != NULL);
   assert(scip != NULL);

   /* get heuristic's data */
   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);

   /* initialize data */
   heurdata->nselectedsols = 0;
   heurdata->ncalls = 0;
   heurdata->nlastsols = 0;
   heurdata->usednodes = 0;
   heurdata->prevlastsol = NULL;
   heurdata->prevbestsol = NULL;
   heurdata->randseed = 0;
   heurdata->lasttuple = NULL;
   heurdata->nfailures = 0;
   heurdata->prevnsols = 0;

   /* initialize hash table */
   SCIP_CALL( SCIPhashtableCreate(&heurdata->hashtable, SCIPblkmem(scip), HASHSIZE_SOLS,
         hashGetKeySols, hashKeyEqSols, hashKeyValSols, NULL) );
   assert(heurdata->hashtable != NULL );

   return SCIP_OKAY;
}

/** deinitialization method of primal heuristic (called before transformed problem is freed) */
static
SCIP_DECL_HEUREXIT(heurExitRec)
{  /*lint --e{715}*/
   SCIP_HEURDATA* heurdata;
   SOLTUPLE* soltuple;

   assert(heur != NULL);
   assert(scip != NULL);

   /* get heuristic data */
   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);
   soltuple = heurdata->lasttuple;

   /* free all soltuples iteratively */
   while( soltuple != NULL )
   {
      SOLTUPLE* tmp;
      tmp = soltuple->prev;
      SCIPfreeBlockMemoryArray(scip,&soltuple->indices,soltuple->size);
      SCIPfreeBlockMemory(scip,&soltuple);
      soltuple = tmp;
   }

   /* free hash table */
   assert(heurdata->hashtable != NULL );
   SCIPhashtableFree(&heurdata->hashtable);

   return SCIP_OKAY;
}


/** execution method of primal heuristic */
static
SCIP_DECL_HEUREXEC(heurExecRec)
{  /*lint --e{715}*/
   SCIP_HEURDATA* heurdata;                  /* primal heuristic data                               */
   SCIP* subscip;                            /* the subproblem created by rec                 */
   SCIP_HASHMAP* varmapfw;                   /* mapping of SCIP variables to sub-SCIP variables */
   SCIP_EVENTHDLR*       eventhdlr;          /* event handler for LP events                     */

   SCIP_VAR** vars;                          /* original problem variables                        */
   SCIP_VAR** subvars;                       /* subproblem variables                              */
   SCIP_SOL** sols;

   SCIP_Real memorylimit;                    /* memory limit for the subproblem                     */
   SCIP_Real timelimit;                      /* time limit for the subproblem                       */
   SCIP_Real cutoff;                         /* objective cutoff for the subproblem                 */
   SCIP_Real upperbound;
   SCIP_Bool success;
   SCIP_Longint nsols;
   SCIP_Longint nstallnodes;                 /* node limit for the subproblem                       */

   int* selection;                           /* pool of solutions rec uses                    */
   int nvars;                                /* number of original problem's variables              */
   int nbinvars;
   int nintvars;
   int nusedsols;
   int i;
   SCIP_RETCODE retcode;

   assert(heur != NULL);
   assert(scip != NULL);
   assert(result != NULL);

   /* get heuristic's data */
   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);
   nusedsols = heurdata->nusedsols;

   *result = SCIP_DELAYED;
   //printf("nsols: %d\n", (int) SCIPgetNSolsFound(scip));
   nsols = SCIPgetNSolsFound(scip);

   /* only call heuristic, if sufficiently many solutions are available */
   if( SCIPgetNSols(scip) < nusedsols  )
      return SCIP_OKAY;

   *result = SCIP_DIDNOTRUN;

   sols = SCIPgetSols(scip);
   assert(sols != NULL);

   /* if a best solution has been found, the heuristic should not be additonally delayed */
   if( sols[0] != heurdata->prevbestsol )
      heurdata->nfailures = 0;


   if( SCIPisLT(scip, nsols, heurdata->nlastsols + heurdata->nwaitingsols + heurdata->nfailures) && heurdata->ncalls > 0 )
   {
      //printf("Recombination suspending\n");
      return SCIP_OKAY;
   }

   heurdata->ncalls++;

   if( SCIPgetNNodes(scip) > 1 )
      SCIPheurSetTimingmask(heur, HEUR_TIMING);
   //printf("rec nfailures:  %d\n", heurdata->nfailures);

   nstallnodes = heurdata->maxnodes;

   if( SCIPisStopped(scip) )
   {
      *result = SCIP_DIDNOTRUN;
      return SCIP_OKAY;
   }

   *result = SCIP_DIDNOTFIND;

   SCIP_CALL( SCIPgetVarsData(scip, &vars, &nvars, &nbinvars, &nintvars, NULL, NULL) );
   assert(nvars > 0);

   /* initializing the subproblem */
   SCIP_CALL( SCIPcreate(&subscip) );

   /* create the variable mapping hash map */
   SCIP_CALL( SCIPhashmapCreate(&varmapfw, SCIPblkmem(subscip), SCIPcalcHashtableSize(5 * nvars)) );
   success = FALSE;

   eventhdlr = NULL;

   if( heurdata->uselprows )
   {
      char probname[SCIP_MAXSTRLEN];

      /* copy all plugins */
      SCIP_CALL( SCIPincludeDefaultPlugins(subscip) );

      /* get name of the original problem and add the string "_recsub" */
      (void) SCIPsnprintf(probname, SCIP_MAXSTRLEN, "%s_recsub", SCIPgetProbName(scip));

      /* create the subproblem */
      SCIP_CALL( SCIPcreateProb(subscip, probname, NULL, NULL, NULL, NULL, NULL, NULL, NULL) );

      /* copy all variables */
      SCIP_CALL( SCIPcopyVars(scip, subscip, varmapfw, NULL, TRUE) );
   }
   else
   {
      SCIP_CALL( SCIPcopy(scip, subscip, varmapfw, NULL, "rec", TRUE, FALSE, TRUE, &success) );

      if( heurdata->copycuts )
      {
         /* copies all active cuts from cutpool of sourcescip to linear constraints in targetscip */
         SCIP_CALL( SCIPcopyCuts(scip, subscip, varmapfw, NULL, TRUE, NULL) );
      }

      /* create event handler for LP events */
      SCIP_CALL( SCIPincludeEventhdlrBasic(subscip, &eventhdlr, EVENTHDLR_NAME, EVENTHDLR_DESC, eventExecRec, NULL) );
      if( eventhdlr == NULL )
      {
         SCIPerrorMessage("event handler for "HEUR_NAME" heuristic not found.\n");
         return SCIP_PLUGINNOTFOUND;
      }
   }

   SCIP_CALL( SCIPallocBufferArray(scip, &subvars, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &selection, nusedsols) );

   for( i = 0; i < nvars; i++ )
      subvars[i] = (SCIP_VAR*) SCIPhashmapGetImage(varmapfw, vars[i]);

   /* free hash map */
   SCIPhashmapFree(&varmapfw);

   success = FALSE;

   SCIP_CALL( SCIPsetIntParam(subscip, "display/verblevel", 0) );
   /* create a new problem, which fixes variables with same value in a certain set of solutions */
   SCIP_CALL( setupSubproblem(scip, subscip, subvars, selection, heurdata, &success) );

   heurdata->prevbestsol = SCIPgetBestSol(scip);
   heurdata->prevlastsol = sols[nsols - 1];
   heurdata->nlastsols   = SCIPgetNSolsFound(scip);
   /* if the creation of sub-SCIP has been aborted (e.g. due to number of fixings), free sub-SCIP and abort */
   if( !success )
   {
      *result = SCIP_DIDNOTRUN;

      /* this run will be counted as a failure since no new solution tuple could be generated or the neighborhood of the
       * solution was not fruitful in the sense that it was too big
       */
      heurdata->nfailures++;

      goto TERMINATE;
   }

   /* do not abort subproblem on CTRL-C */
   SCIP_CALL( SCIPsetBoolParam(subscip, "misc/catchctrlc", FALSE) );

   /* disable output to console */

#ifdef SCIP_DEBUG
   /* for debugging DINS, enable MIP output */
   SCIP_CALL( SCIPsetIntParam(subscip, "display/verblevel", 5) );
   SCIP_CALL( SCIPsetIntParam(subscip, "display/freq", 1) );
#endif

   /* check whether there is enough time and memory left */
   SCIP_CALL( SCIPgetRealParam(scip, "limits/time", &timelimit) );
   if( !SCIPisInfinity(scip, timelimit) )
      timelimit -= SCIPgetSolvingTime(scip);
   SCIP_CALL( SCIPgetRealParam(scip, "limits/memory", &memorylimit) );

   /* substract the memory already used by the main SCIP and the estimated memory usage of external software */
   if( !SCIPisInfinity(scip, memorylimit) )
   {
      memorylimit -= SCIPgetMemUsed(scip)/1048576.0;
      memorylimit -= SCIPgetMemExternEstim(scip)/1048576.0;
   }

   /* abort if no time is left or not enough memory to create a copy of SCIP, including external memory usage */
   if( timelimit <= 0.0 || memorylimit <= 2.0*SCIPgetMemExternEstim(scip)/1048576.0 )
      goto TERMINATE;

   /* set limits for the subproblem */
   heurdata->nodelimit = nstallnodes;
   SCIP_CALL( SCIPsetLongintParam(subscip, "limits/nodes", nstallnodes) );
   SCIP_CALL( SCIPsetRealParam(subscip, "limits/time", timelimit) );
   SCIP_CALL( SCIPsetRealParam(subscip, "limits/memory", memorylimit) );

   /* deactivate separators and heuristics that use auxiliary SCIP instances */
   SCIP_CALL( SCIPsetSubscipsOff(subscip, TRUE) );
#if 0
   /* disable cutting plane separation */
   SCIP_CALL( SCIPsetSeparating(subscip, SCIP_PARAMSETTING_OFF, TRUE) );
#endif
   /* disable expensive presolving */
   SCIP_CALL( SCIPsetPresolving(subscip, SCIP_PARAMSETTING_FAST, TRUE) );

   /* use best estimate node selection */
   if( SCIPfindNodesel(subscip, "estimate") != NULL && !SCIPisParamFixed(subscip, "nodeselection/estimate/stdpriority") )
   {
      SCIP_CALL( SCIPsetIntParam(subscip, "nodeselection/estimate/stdpriority", INT_MAX/4) );
   }

   /* use inference branching */
   if( SCIPfindBranchrule(subscip, "inference") != NULL && !SCIPisParamFixed(subscip, "branching/inference/priority") )
   {
      SCIP_CALL( SCIPsetIntParam(subscip, "branching/inference/priority", INT_MAX/4) );
   }

   /* disable conflict analysis */
   if( !SCIPisParamFixed(subscip, "conflict/useprop") )
   {
      SCIP_CALL( SCIPsetBoolParam(subscip, "conflict/useprop", FALSE) );
   }
   if( !SCIPisParamFixed(subscip, "conflict/useinflp") )
   {
      SCIP_CALL( SCIPsetBoolParam(subscip, "conflict/useinflp", FALSE) );
   }
   if( !SCIPisParamFixed(subscip, "conflict/useboundlp") )
   {
      SCIP_CALL( SCIPsetBoolParam(subscip, "conflict/useboundlp", FALSE) );
   }
   if( !SCIPisParamFixed(subscip, "conflict/usesb") )
   {
      SCIP_CALL( SCIPsetBoolParam(subscip, "conflict/usesb", FALSE) );
   }
   if( !SCIPisParamFixed(subscip, "conflict/usepseudo") )
   {
      SCIP_CALL( SCIPsetBoolParam(subscip, "conflict/usepseudo", FALSE) );
   }

   /* employ a limit on the number of enforcement rounds in the quadratic constraint handler; this fixes the issue that
    * sometimes the quadratic constraint handler needs hundreds or thousands of enforcement rounds to determine the
    * feasibility status of a single node without fractional branching candidates by separation (namely for uflquad
    * instances); however, the solution status of the sub-SCIP might get corrupted by this; hence no deductions shall be
    * made for the original SCIP
    */
   if( SCIPfindConshdlr(subscip, "quadratic") != NULL && !SCIPisParamFixed(subscip, "constraints/quadratic/enfolplimit") )
   {
      SCIP_CALL( SCIPsetIntParam(subscip, "constraints/quadratic/enfolplimit", 500) );
   }

   /* add an objective cutoff */
   cutoff = SCIPinfinity(scip);
   assert(!SCIPisInfinity(scip, SCIPgetUpperbound(scip)));

   upperbound = SCIPgetUpperbound(scip) - SCIPsumepsilon(scip);
   if( !SCIPisInfinity(scip, -1.0 * SCIPgetLowerbound(scip)) )
   {
      cutoff = (1-heurdata->minimprove)*SCIPgetUpperbound(scip) + heurdata->minimprove*SCIPgetLowerbound(scip);
   }
   else
   {
      if( SCIPgetUpperbound ( scip ) >= 0 )
         cutoff = ( 1 - heurdata->minimprove ) * SCIPgetUpperbound ( scip );
      else
         cutoff = ( 1 + heurdata->minimprove ) * SCIPgetUpperbound ( scip );
   }
   cutoff = MIN(upperbound, cutoff );
   SCIP_CALL( SCIPsetObjlimit(subscip, cutoff) );

   /* permute the subproblem to increase diversification */
   if( heurdata->permute )
      SCIP_CALL( SCIPpermuteProb(subscip, (unsigned int) SCIPheurGetNCalls(heur), TRUE, TRUE, TRUE, TRUE, TRUE) );

   /* catch LP events of sub-SCIP */
   if( !heurdata->uselprows )
   {
      assert(eventhdlr != NULL);

      SCIP_CALL( SCIPtransformProb(subscip) );
      SCIP_CALL( SCIPcatchEvent(subscip, SCIP_EVENTTYPE_LPSOLVED, eventhdlr, (SCIP_EVENTDATA*) heurdata, NULL) );
   }

   /* solve the subproblem */
   SCIPdebugMessage("Solve Rec subMIP\n");
   retcode = SCIPsolve(subscip);

   /* drop LP events of sub-SCIP */
   if( !heurdata->uselprows )
   {
      assert(eventhdlr != NULL);

      SCIP_CALL( SCIPdropEvent(subscip, SCIP_EVENTTYPE_LPSOLVED, eventhdlr, (SCIP_EVENTDATA*) heurdata, -1) );
   }

   /* Errors in solving the subproblem should not kill the overall solving process.
    * Hence, the return code is caught and a warning is printed, only in debug mode, SCIP will stop. */
   if( retcode != SCIP_OKAY )
   {
#ifndef NDEBUG
      SCIP_CALL( retcode );
#endif
      SCIPwarningMessage(scip, "Error while solving subproblem in Rec heuristic; sub-SCIP terminated with code <%d>\n", retcode);
   }

   /* print solving statistics of subproblem if we are in SCIP's debug mode */
   SCIPdebug( SCIP_CALL( SCIPprintStatistics(subscip, NULL) ) );

   heurdata->usednodes += SCIPgetNNodes(subscip);

   /* check, whether a solution was found */
   if( SCIPgetNSols(subscip) > 0 )
   {
      SCIP_SOL** subsols;
      int nsubsols;
      int solindex;                             /* index of the solution created by rec          */

      /* check, whether a solution was found;
       * due to numerics, it might happen that not all solutions are feasible -> try all solutions until one was accepted */
      nsubsols = SCIPgetNSols(subscip);
      subsols = SCIPgetSols(subscip);
      success = FALSE;
      solindex = -1;
      for( i = 0; i < nsubsols && !success; ++i )
      {
         SCIP_CALL( createNewSol(scip, subscip, subvars, heur, subsols[i], &solindex, &success) );
      }

      if( success )
      {
         int tmp;

         assert(solindex != -1);

         *result = SCIP_FOUNDSOL;

         /* insert all crossings of the new solution and (nusedsols-1) of its parents into the hashtable
          * in order to avoid incest ;)
          */
         for( i = 0; i < heurdata->nselectedsols; i++ )
         {
            SOLTUPLE* elem;
            tmp = selection[i];
            selection[i] = solindex;

            SCIP_CALL( createSolTuple(scip, &elem, selection, heurdata->nselectedsols, heurdata) );
            SCIP_CALL( SCIPhashtableInsert(heurdata->hashtable, elem) );
            selection[i] = tmp;
         }

         /* if solution was among the best ones, rec should not be called until another good solution was found */
         heurdata->prevbestsol = SCIPgetBestSol(scip);
         heurdata->prevlastsol = SCIPgetSols(scip)[nsols - 1];
         heurdata->nlastsols   = SCIPgetNSolsFound(scip);

      }

      /* if solution is not better then incumbent or could not be added to problem => run is counted as a failure */
      if( !success || solindex != SCIPsolGetIndex(SCIPgetBestSol(scip)) )
         heurdata->nfailures++;
   }
   else
   {
      /* if no new solution has been found, run is considered a failure */
      heurdata->nfailures++;
   }

 TERMINATE:
   /* free subproblem */
   SCIPfreeBufferArray(scip, &selection);
   SCIPfreeBufferArray(scip, &subvars);
   SCIP_CALL( SCIPfree(&subscip) );

   return SCIP_OKAY;
}


/*
 * primal heuristic specific interface methods
 */

/** creates the rec primal heuristic and includes it in SCIP */
SCIP_RETCODE SCIPincludeHeurRec(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_HEURDATA* heurdata;
   SCIP_HEUR* heur;

   /* create Rec primal heuristic data */
   SCIP_CALL( SCIPallocMemory(scip, &heurdata) );

   /* include primal heuristic */
   SCIP_CALL( SCIPincludeHeurBasic(scip, &heur,
         HEUR_NAME, HEUR_DESC, HEUR_DISPCHAR, HEUR_PRIORITY, HEUR_FREQ, HEUR_FREQOFS,
         HEUR_MAXDEPTH, HEUR_TIMING, HEUR_USESSUBSCIP, heurExecRec, heurdata) );

   assert(heur != NULL);

   /* set non-NULL pointers to callback methods */
   SCIP_CALL( SCIPsetHeurCopy(scip, heur, heurCopyRec) );
   SCIP_CALL( SCIPsetHeurFree(scip, heur, heurFreeRec) );
   SCIP_CALL( SCIPsetHeurInit(scip, heur, heurInitRec) );
   SCIP_CALL( SCIPsetHeurExit(scip, heur, heurExitRec) );

   /* add rec primal heuristic parameters */

   SCIP_CALL( SCIPaddLongintParam(scip, "heuristics/"HEUR_NAME"/nodesofs",
         "number of nodes added to the contingent of the total nodes",
         &heurdata->nodesofs, FALSE, DEFAULT_NODESOFS, 0LL, SCIP_LONGINT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddLongintParam(scip, "heuristics/"HEUR_NAME"/maxnodes",
         "maximum number of nodes to regard in the subproblem",
         &heurdata->maxnodes, TRUE, DEFAULT_MAXNODES, 0LL, SCIP_LONGINT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddLongintParam(scip, "heuristics/"HEUR_NAME"/minnodes",
         "minimum number of nodes required to start the subproblem",
         &heurdata->minnodes, TRUE, DEFAULT_MINNODES, 0LL, SCIP_LONGINT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(scip, "heuristics/"HEUR_NAME"/nusedsols",
         "number of solutions to be taken into account",
         &heurdata->nusedsols, FALSE, DEFAULT_NUSEDSOLS, 2, INT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(scip, "heuristics/"HEUR_NAME"/nwaitingsols",
         "number of solutions to be in abeyance",
         &heurdata->nwaitingsols, FALSE, DEFAULT_NWAITINGSOLS, 1, INT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddLongintParam(scip, "heuristics/"HEUR_NAME"/nwaitingnodes",
         "number of nodes without incumbent change that heuristic should wait",
         &heurdata->nwaitingnodes, TRUE, DEFAULT_NWAITINGNODES, 0LL, SCIP_LONGINT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip, "heuristics/"HEUR_NAME"/nodesquot",
         "contingent of sub problem nodes in relation to the number of nodes of the original problem",
         &heurdata->nodesquot, FALSE, DEFAULT_NODESQUOT, 0.0, 1.0, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip, "heuristics/"HEUR_NAME"/minfixingrate",
         "minimum percentage of integer variables that have to be fixed",
         &heurdata->minfixingrate, FALSE, DEFAULT_MINFIXINGRATE, 0.0, 1.0, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip, "heuristics/"HEUR_NAME"/minimprove",
         "factor by which Rec should at least improve the incumbent",
         &heurdata->minimprove, TRUE, DEFAULT_MINIMPROVE, 0.0, 1.0, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip, "heuristics/"HEUR_NAME"/lplimfac",
         "factor by which the limit on the number of LP depends on the node limit",
         &heurdata->lplimfac, TRUE, DEFAULT_LPLIMFAC, 1.0, SCIP_REAL_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "heuristics/"HEUR_NAME"/dontwaitatroot",
         "should the nwaitingnodes parameter be ignored at the root node?",
         &heurdata->dontwaitatroot, TRUE, DEFAULT_DONTWAITATROOT, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "heuristics/"HEUR_NAME"/uselprows",
         "should subproblem be created out of the rows in the LP rows?",
         &heurdata->uselprows, TRUE, DEFAULT_USELPROWS, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "heuristics/"HEUR_NAME"/copycuts",
         "if uselprows == FALSE, should all active cuts from cutpool be copied to constraints in subproblem?",
         &heurdata->copycuts, TRUE, DEFAULT_COPYCUTS, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "heuristics/"HEUR_NAME"/permute",
         "should the subproblem be permuted to increase diversification?",
         &heurdata->permute, TRUE, DEFAULT_PERMUTE, NULL, NULL) );

   return SCIP_OKAY;
}
