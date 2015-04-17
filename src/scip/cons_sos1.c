/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2015 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   cons_sos1.c
 * @brief  constraint handler for SOS type 1 constraints
 * @author Tobias Fischer
 * @author Marc Pfetsch
 *
 * A specially ordered set of type 1 (SOS1) is a sequence of variables such that at most one
 * variable is nonzero. The special case of two variables arises, for instance, from equilibrium or
 * complementary conditions like \f$x \cdot y = 0\f$. Note that it is in principle allowed that a
 * variables appears twice, but it then can be fixed to 0.
 *
 * This implementation of this constraint handler is based on classical ideas, see e.g.@n
 *  "Special Facilities in General Mathematical Programming System for
 *  Non-Convex Problems Using Ordered Sets of Variables"@n
 *  E. Beale and J. Tomlin, Proc. 5th IFORS Conference, 447-454 (1970)
 *
 *
 * The order of the variables is determined as follows:
 *
 * - If the constraint is created with SCIPcreateConsSOS1() and weights are given, the weights
 *   determine the order (decreasing weights). Additional variables can be added with
 *   SCIPaddVarSOS1(), which adds a variable with given weight.
 *
 * - If an empty constraint is created and then variables are added with SCIPaddVarSOS1(), weights
 *   are needed and stored.
 *
 * - All other calls ignore the weights, i.e., if a nonempty constraint is created or variables are
 *   added with SCIPappendVarSOS1().
 *
 * The validity of the SOS1 constraint can be enforced by different branching rules:
 *
 * - If classical SOS branching is used, branching is performed on only one SOS1 constraint. Depending on the parameters,
 *   there are two ways to choose this branching constraint. Either the constraint with the most number of nonzeros
 *   or the one with the largest nonzero-variable weight. The later version allows the user to specify
 *   an order for the branching importance of the constraints. Constraint branching can also be turned off.
 *
 * - Another way is to branch on the neighborhood of a single variable @p i, i.e., in one branch \f$x_i\f$ is fixed to zero
 *   and in the other its neighbors.
 *
 * - If bipartite branching is used, then we branch using complete bipartite subgraphs.
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>

#include "scip/cons_sos1.h"
#include "scip/cons_linear.h"
#include "scip/cons_setppc.h"
#include "scip/pub_misc.h"
#include "scip/misc.h"
#include "scip/struct_misc.h"
#include "tclique/tclique.h"
#include <string.h>
#include <ctype.h>


/* constraint handler properties */
#define CONSHDLR_NAME          "SOS1"
#define CONSHDLR_DESC          "SOS1 constraint handler"
#define CONSHDLR_SEPAPRIORITY   -900000 /**< priority of the constraint handler for separation */
#define CONSHDLR_ENFOPRIORITY       100 /**< priority of the constraint handler for constraint enforcing */
#define CONSHDLR_CHECKPRIORITY      -10 /**< priority of the constraint handler for checking feasibility */
#define CONSHDLR_SEPAFREQ            10 /**< frequency for separating cuts; zero means to separate only in the root node */
#define CONSHDLR_PROPFREQ             1 /**< frequency for propagating domains; zero means only preprocessing propagation */
#define CONSHDLR_EAGERFREQ          100 /**< frequency for using all instead of only the useful constraints in separation,
                                         *   propagation and enforcement, -1 for no eager evaluations, 0 for first only */
#define CONSHDLR_MAXPREROUNDS        -1 /**< maximal number of presolving rounds the constraint handler participates in (-1: no limit) */
#define CONSHDLR_DELAYSEPA        FALSE /**< should separation method be delayed, if other separators found cuts? */
#define CONSHDLR_DELAYPROP        FALSE /**< should propagation method be delayed, if other propagators found reductions? */
#define CONSHDLR_DELAYPRESOL       TRUE /**< should presolving method be delayed, if other presolvers found reductions? */
#define CONSHDLR_NEEDSCONS         TRUE /**< should the constraint handler be skipped, if no constraints are available? */

/* presolving */
#define DEFAULT_MAXEXTENSIONS         1 /**< maximal number of extensions that will be computed for each SOS1 constraint */
#define DEFAULT_MAXCONSDELAYEXT      -1 /**< delay clique extension if number of SOS1 constraints is larger than predefined value (-1: no limit) */
#define DEFAULT_MAXTIGHTENBDS        -1 /**< maximal number of bound tightening rounds per presolving round (-1: no limit) */
#define DEFAULT_UPDATECONFLPRESOL FALSE /**< if TRUE then update conflict graph during presolving procedure */

/* propagation */
#define DEFAULT_CONFLICTPROP      TRUE /**< whether to use conflict graph propagation */
#define DEFAULT_SOSCONSPROP      FALSE /**< whether to use SOS1 constraint propagation */

/* separation */
#define DEFAULT_SEPAFROMSOS1      FALSE /**< if TRUE separate bound inequalities from initial SOS1 constraints */
#define DEFAULT_SEPAFROMGRAPH      TRUE /**< if TRUE separate bound inequalities from the conflict graph */
#define DEFAULT_BOUNDCUTSDEPTH       40 /**< node depth of separating bound cuts (-1: no limit) */
#define DEFAULT_MAXBOUNDCUTS         50 /**< maximal number of bound cuts separated per branching node */
#define DEFAULT_MAXBOUNDCUTSROOT    150 /**< maximal number of bound cuts separated per iteration in the root node */
#define DEFAULT_STRTHENBOUNDCUTS   TRUE /**< if TRUE then bound cuts are strengthened in case bound variables are available */

#define CONSHDLR_PROP_TIMING       SCIP_PROPTIMING_BEFORELP

/* event handler properties */
#define EVENTHDLR_NAME         "SOS1"
#define EVENTHDLR_DESC         "bound change event handler for SOS1 constraints"


/** constraint data for SOS1 constraints */
struct SCIP_ConsData
{
   int                   nvars;              /**< number of variables in the constraint */
   int                   maxvars;            /**< maximal number of variables (= size of storage) */
   int                   nfixednonzeros;     /**< number of variables fixed to be nonzero */
   SCIP_Bool             local;              /**< TRUE if constraint is only valid locally */
   SCIP_VAR**            vars;               /**< variables in constraint */
   SCIP_ROW*             rowlb;              /**< row corresponding to lower bounds, or NULL if not yet created */
   SCIP_ROW*             rowub;              /**< row corresponding to upper bounds, or NULL if not yet created */
   SCIP_Real*            weights;            /**< weights determining the order (ascending), or NULL if not used */
};


/** node data of a given node in the conflict graph */
struct SCIP_NodeData
{
   SCIP_VAR*             var;                /**< variable belonging to node */
   SCIP_VAR*             lbboundvar;         /**< bound variable @p z from constraint \f$x \geq \mu \cdot z\f$ (or NULL if not existent) */
   SCIP_VAR*             ubboundvar;         /**< bound variable @p z from constraint \f$x \leq \mu \cdot z\f$ (or NULL if not existent) */
   SCIP_Real             lbboundcoef;        /**< value \f$\mu\f$ from constraint \f$x \geq \mu z \f$ (0.0 if not existent) */
   SCIP_Real             ubboundcoef;        /**< value \f$\mu\f$ from constraint \f$x \leq \mu z \f$ (0.0 if not existent) */
   SCIP_Bool             lbboundcomp;        /**< TRUE if the nodes from the connected component of the conflict graph the given node belongs to
                                              *   all have the same lower bound variable */
   SCIP_Bool             ubboundcomp;        /**< TRUE if the nodes from the connected component of the conflict graph the given node belongs to
                                              *   all have the same lower bound variable */
};
typedef struct SCIP_NodeData SCIP_NODEDATA;


/** successor data of a given nodes successor in the implication graph */
struct SCIP_SuccData
{
   SCIP_Real             lbimpl;             /**< lower bound implication */
   SCIP_Real             ubimpl;             /**< upper bound implication */
};
typedef struct SCIP_SuccData SCIP_SUCCDATA;


/** tclique data for bound cut generation */
struct TCLIQUE_Data
{
   SCIP*                 scip;               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr;           /**< SOS1 constraint handler */
   SCIP_DIGRAPH*         conflictgraph;      /**< conflict graph */
   SCIP_SOL*             sol;                /**< LP solution to be separated (or NULL) */
   SCIP_Real             scaleval;           /**< factor for scaling weights */
   int                   ncuts;              /**< number of bound cuts found in this iteration */
   int                   nboundcuts;         /**< number of bound cuts found so far */
   int                   maxboundcuts;       /**< maximal number of clique cuts separated per separation round (-1: no limit) */
   SCIP_Bool             strthenboundcuts;   /**< if TRUE then bound cuts are strengthened in case bound variables are available */
};


/** SOS1 constraint handler data */
struct SCIP_ConshdlrData
{
   /* conflict graph */
   SCIP_DIGRAPH*         conflictgraph;      /**< conflict graph */
   SCIP_DIGRAPH*         localconflicts;     /**< local conflicts */
   SCIP_Bool             isconflocal;        /**< if TRUE then local conflicts are present and conflict graph has to be updated for each node */
   SCIP_HASHMAP*         varhash;            /**< hash map from variable to node in the conflict graph */
   int                   nsos1vars;          /**< number of problem variables that are involved in at least one SOS1 constraint */
   /* presolving */
   int                   cntextsos1;         /**< counts number of extended SOS1 constraints */
   int                   maxextensions;      /**< maximal number of extensions that will be computed for each SOS1 constraint */
   int                   maxconsdelayext;    /**< delay clique extension if number of sos1 constraints is larger than predefined value (-1: no limit) */
   int                   maxtightenbds;      /**< maximal number of bound tightening rounds per presolving round (-1: no limit) */
   SCIP_Bool             updateconflpresol;  /**< if TRUE then update conflict graph during presolving procedure */
   /* propagation */
   SCIP_Bool             conflictprop;       /**< whether to use conflict graph propagation */
   SCIP_Bool             sosconsprop;        /**< whether to use SOS1 constraint propagation */
   /* branching */
   SCIP_Bool             branchsos;          /**< Branch on SOS condition in enforcing? */
   SCIP_Bool             branchnonzeros;     /**< Branch on SOS cons. with most number of nonzeros? */
   SCIP_Bool             branchweight;       /**< Branch on SOS cons. with highest nonzero-variable weight for branching - needs branchnonzeros to be false */
   /* separation */
   SCIP_Bool             sepafromsos1;       /**< if TRUE separate bound inequalities from initial SOS1 constraints */
   SCIP_Bool             sepafromgraph;      /**< if TRUE separate bound inequalities from the conflict graph */
   TCLIQUE_GRAPH*        tcliquegraph;       /**< tclique graph data structure */
   TCLIQUE_DATA*         tcliquedata;        /**< tclique data */
   int                   boundcutsdepth;     /**< node depth of separating bound cuts (-1: no limit) */
   int                   maxboundcuts;       /**< maximal number of bound cuts separated per branching node */
   int                   maxboundcutsroot;   /**< maximal number of bound cuts separated per iteration in the root node */
   int                   nboundcuts;         /**< number of bound cuts found so far */
   SCIP_Bool             strthenboundcuts;   /**< if TRUE then bound cuts are strengthened in case bound variables are available */
   /* event handler */
   SCIP_EVENTHDLR*       eventhdlr;          /**< event handler for bound change events */
};


/** returns whether two vertices are adjacent (connected) */
static
SCIP_Bool isConnectedSOS1(
   SCIP_Bool**           adjacencymatrix,    /**< adjacency matrix of conflict graph (lower half) (or NULL if an adjacencymatrix is not at hand) */
   SCIP_DIGRAPH*         conflictgraph,      /**< conflict graph (or NULL if an adjacencymatrix is at hand) */
   int                   vertex1,            /**< first vertex */
   int                   vertex2             /**< second vertex */
   )
{
   assert( adjacencymatrix != NULL || conflictgraph != NULL );

   /* we do not allow self-loops */
   if ( vertex1 == vertex2 )
      return FALSE;

   /* for debugging */
   if ( adjacencymatrix == NULL )
   {
      int succvertex;
      int* succ;
      int nsucc1;
      int nsucc2;
      int j;

      nsucc1 = SCIPdigraphGetNSuccessors(conflictgraph, vertex1);
      nsucc2 = SCIPdigraphGetNSuccessors(conflictgraph, vertex2);

      if ( nsucc1 < 1 || nsucc2 < 1 )
         return FALSE;

      if ( nsucc1 > nsucc2 )
      {
         SCIPswapInts(&vertex1, &vertex2);
         SCIPswapInts(&nsucc1, &nsucc2);
      }

      succ = SCIPdigraphGetSuccessors(conflictgraph, vertex1);
      SCIPsortInt(succ, nsucc1);

      for (j = 0; j < nsucc1; ++j)
      {
         succvertex = succ[j];
         if ( succvertex == vertex2 )
            return TRUE;
         else if ( succvertex > vertex2 )
            return FALSE;
      }
   }
   else
   {
      if ( vertex1 < vertex2 )
         return adjacencymatrix[vertex2][vertex1];
      else
         return adjacencymatrix[vertex1][vertex2];
   }

   return FALSE;
}


/** fix variable in given node to 0 or add constraint if variable is multi-aggregated */
static
SCIP_RETCODE fixVariableZeroNode(
   SCIP*                 scip,               /**< SCIP pointer */
   SCIP_VAR*             var,                /**< variable to be fixed to 0*/
   SCIP_NODE*            node,               /**< node */
   SCIP_Bool*            infeasible          /**< if fixing is infeasible */
   )
{
   /* if variable cannot be nonzero */
   *infeasible = FALSE;
   if ( SCIPisFeasPositive(scip, SCIPvarGetLbLocal(var)) || SCIPisFeasNegative(scip, SCIPvarGetUbLocal(var)) )
   {
      *infeasible = TRUE;
      return SCIP_OKAY;
   }

   /* if variable is multi-aggregated */
   if ( SCIPvarGetStatus(var) == SCIP_VARSTATUS_MULTAGGR )
   {
      SCIP_CONS* cons;
      SCIP_Real val;

      val = 1.0;

      if ( ! SCIPisFeasZero(scip, SCIPvarGetLbLocal(var)) || ! SCIPisFeasZero(scip, SCIPvarGetUbLocal(var)) )
      {
         SCIPdebugMessage("creating constraint to force multi-aggregated variable <%s> to 0.\n", SCIPvarGetName(var));
         /* we have to insert a local constraint var = 0 */
         SCIP_CALL( SCIPcreateConsLinear(scip, &cons, "branch", 1, &var, &val, 0.0, 0.0, TRUE, TRUE, TRUE, TRUE, TRUE,
               TRUE, FALSE, FALSE, FALSE, FALSE) );
         SCIP_CALL( SCIPaddConsNode(scip, node, cons, NULL) );
         SCIP_CALL( SCIPreleaseCons(scip, &cons) );
      }
   }
   else
   {
      if ( ! SCIPisFeasZero(scip, SCIPvarGetLbLocal(var)) )
         SCIP_CALL( SCIPchgVarLbNode(scip, node, var, 0.0) );
      if ( ! SCIPisFeasZero(scip, SCIPvarGetUbLocal(var)) )
         SCIP_CALL( SCIPchgVarUbNode(scip, node, var, 0.0) );
   }

   return SCIP_OKAY;
}


/** fix variable in local node to 0, and return whether the operation was feasible
 *
 *  @note We do not add a linear constraint if the variable is multi-aggregated as in
 *  fixVariableZeroNode(), since this would be too time consuming.
 */
static
SCIP_RETCODE inferVariableZero(
   SCIP*                 scip,               /**< SCIP pointer */
   SCIP_VAR*             var,                /**< variable to be fixed to 0*/
   SCIP_CONS*            cons,               /**< constraint */
   int                   inferinfo,          /**< info for reverse prop. */
   SCIP_Bool*            infeasible,         /**< if fixing is infeasible */
   SCIP_Bool*            tightened,          /**< if fixing was performed */
   SCIP_Bool*            success             /**< whether fixing was successful, i.e., variable is not multi-aggregated */
   )
{
   *infeasible = FALSE;
   *tightened = FALSE;
   *success = FALSE;

   /* if variable cannot be nonzero */
   if ( SCIPisFeasPositive(scip, SCIPvarGetLbLocal(var)) || SCIPisFeasNegative(scip, SCIPvarGetUbLocal(var)) )
   {
      *infeasible = TRUE;
      return SCIP_OKAY;
   }

   /* directly fix variable if it is not multi-aggregated */
   if ( SCIPvarGetStatus(var) != SCIP_VARSTATUS_MULTAGGR )
   {
      SCIP_Bool tighten;

      /* fix lower bound */
      SCIP_CALL( SCIPinferVarLbCons(scip, var, 0.0, cons, inferinfo, FALSE, infeasible, &tighten) );
      *tightened = *tightened || tighten;

      /* fix upper bound */
      SCIP_CALL( SCIPinferVarUbCons(scip, var, 0.0, cons, inferinfo, FALSE, infeasible, &tighten) );
      *tightened = *tightened || tighten;

      *success = TRUE;
   }

   return SCIP_OKAY;
}


/** add lock on variable */
static
SCIP_RETCODE lockVariableSOS1(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint */
   SCIP_VAR*             var                 /**< variable */
   )
{
   assert( scip != NULL );
   assert( cons != NULL );
   assert( var != NULL );

   /* rounding down == bad if lb < 0, rounding up == bad if ub > 0 */
   SCIP_CALL( SCIPlockVarCons(scip, var, cons, SCIPisFeasNegative(scip, SCIPvarGetLbLocal(var)), SCIPisFeasPositive(scip, SCIPvarGetUbLocal(var))) );

   return SCIP_OKAY;
}


/* remove lock on variable */
static
SCIP_RETCODE unlockVariableSOS1(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint */
   SCIP_VAR*             var                 /**< variable */
   )
{
   assert( scip != NULL );
   assert( cons != NULL );
   assert( var != NULL );

   /* rounding down == bad if lb < 0, rounding up == bad if ub > 0 */
   SCIP_CALL( SCIPunlockVarCons(scip, var, cons, SCIPisFeasNegative(scip, SCIPvarGetLbLocal(var)), SCIPisFeasPositive(scip, SCIPvarGetUbLocal(var))) );

   return SCIP_OKAY;
}


/** ensures that the vars and weights array can store at least num entries */
static
SCIP_RETCODE consdataEnsurevarsSizeSOS1(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSDATA*        consdata,           /**< constraint data */
   int                   num,                /**< minimum number of entries to store */
   SCIP_Bool             reserveWeights      /**< whether the weights array is handled */
   )
{
   assert( consdata != NULL );
   assert( consdata->nvars <= consdata->maxvars );

   if ( num > consdata->maxvars )
   {
      int newsize;

      newsize = SCIPcalcMemGrowSize(scip, num);
      SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &consdata->vars, consdata->maxvars, newsize) );
      if ( reserveWeights )
         SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &consdata->weights, consdata->maxvars, newsize) );
      consdata->maxvars = newsize;
   }
   assert( num <= consdata->maxvars );

   return SCIP_OKAY;
}


/** handle new variable */
static
SCIP_RETCODE handleNewVariableSOS1(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint */
   SCIP_CONSDATA*        consdata,           /**< constraint data */
   SCIP_VAR*             var,                /**< variable */
   SCIP_Bool             transformed         /**< whether original variable was transformed */
   )
{
   assert( scip != NULL );
   assert( cons != NULL );
   assert( consdata != NULL );
   assert( var != NULL );

   /* if we are in transformed problem, catch the variable's events */
   if ( transformed )
   {
      SCIP_CONSHDLR* conshdlr;
      SCIP_CONSHDLRDATA* conshdlrdata;

      /* get event handler */
      conshdlr = SCIPconsGetHdlr(cons);
      conshdlrdata = SCIPconshdlrGetData(conshdlr);
      assert( conshdlrdata != NULL );
      assert( conshdlrdata->eventhdlr != NULL );

      /* catch bound change events of variable */
      SCIP_CALL( SCIPcatchVarEvent(scip, var, SCIP_EVENTTYPE_BOUNDCHANGED, conshdlrdata->eventhdlr,
            (SCIP_EVENTDATA*)consdata, NULL) );

      /* if the variable if fixed to nonzero */
      assert( consdata->nfixednonzeros >= 0 );
      if ( SCIPisFeasPositive(scip, SCIPvarGetLbLocal(var)) || SCIPisFeasNegative(scip, SCIPvarGetUbLocal(var)) )
         ++consdata->nfixednonzeros;
   }

   /* install the rounding locks for the new variable */
   SCIP_CALL( lockVariableSOS1(scip, cons, var) );

   /* branching on multiaggregated variables does not seem to work well, so avoid it */
   SCIP_CALL( SCIPmarkDoNotMultaggrVar(scip, var) );

   /* add the new coefficient to the upper bound LP row, if necessary */
   if ( consdata->rowub != NULL && ! SCIPisInfinity(scip, SCIPvarGetUbGlobal(var)) && ! SCIPisZero(scip, SCIPvarGetUbGlobal(var)) )
   {
      SCIP_CALL( SCIPaddVarToRow(scip, consdata->rowub, var, 1.0/SCIPvarGetUbGlobal(var)) );
   }

   /* add the new coefficient to the lower bound LP row, if necessary */
   if ( consdata->rowlb != NULL && ! SCIPisInfinity(scip, SCIPvarGetLbGlobal(var)) && ! SCIPisZero(scip, SCIPvarGetLbGlobal(var)) )
   {
      SCIP_CALL( SCIPaddVarToRow(scip, consdata->rowlb, var, 1.0/SCIPvarGetLbGlobal(var)) );
   }

   return SCIP_OKAY;
}


/** adds a variable to an SOS1 constraint, at position given by weight - ascending order */
static
SCIP_RETCODE addVarSOS1(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint */
   SCIP_VAR*             var,                /**< variable to add to the constraint */
   SCIP_Real             weight              /**< weight to determine position */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_Bool transformed;
   int pos;
   int j;

   assert( var != NULL );
   assert( cons != NULL );

   consdata = SCIPconsGetData(cons);
   assert( consdata != NULL );

   if ( consdata->weights == NULL && consdata->maxvars > 0 )
   {
      SCIPerrorMessage("cannot add variable to SOS1 constraint <%s> that does not contain weights.\n", SCIPconsGetName(cons));
      return SCIP_INVALIDCALL;
   }

   /* are we in the transformed problem? */
   transformed = SCIPconsIsTransformed(cons);

   /* always use transformed variables in transformed constraints */
   if ( transformed )
   {
      SCIP_CALL( SCIPgetTransformedVar(scip, var, &var) );
   }
   assert( var != NULL );
   assert( transformed == SCIPvarIsTransformed(var) );

   SCIP_CALL( consdataEnsurevarsSizeSOS1(scip, consdata, consdata->nvars + 1, TRUE) );
   assert( consdata->weights != NULL );
   assert( consdata->maxvars >= consdata->nvars+1 );

   /* find variable position */
   for (pos = 0; pos < consdata->nvars; ++pos)
   {
      if ( consdata->weights[pos] > weight )
         break;
   }
   assert( 0 <= pos && pos <= consdata->nvars );

   /* move other variables, if necessary */
   for (j = consdata->nvars; j > pos; --j)
   {
      consdata->vars[j] = consdata->vars[j-1];
      consdata->weights[j] = consdata->weights[j-1];
   }

   /* insert variable */
   consdata->vars[pos] = var;
   consdata->weights[pos] = weight;
   ++consdata->nvars;

   /* handle the new variable */
   SCIP_CALL( handleNewVariableSOS1(scip, cons, consdata, var, transformed) );

   return SCIP_OKAY;
}


/** appends a variable to an SOS1 constraint */
static
SCIP_RETCODE appendVarSOS1(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint */
   SCIP_VAR*             var                 /**< variable to add to the constraint */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_Bool transformed;

   assert( var != NULL );
   assert( cons != NULL );

   consdata = SCIPconsGetData(cons);
   assert( consdata != NULL );

   /* are we in the transformed problem? */
   transformed = SCIPconsIsTransformed(cons);

   /* always use transformed variables in transformed constraints */
   if ( transformed )
   {
      SCIP_CALL( SCIPgetTransformedVar(scip, var, &var) );
   }
   assert( var != NULL );
   assert( transformed == SCIPvarIsTransformed(var) );

   SCIP_CALL( consdataEnsurevarsSizeSOS1(scip, consdata, consdata->nvars + 1, FALSE) );

   /* insert variable */
   consdata->vars[consdata->nvars] = var;
   assert( consdata->weights != NULL || consdata->nvars > 0 );
   if ( consdata->weights != NULL && consdata->nvars > 0 )
      consdata->weights[consdata->nvars] = consdata->weights[consdata->nvars-1] + 1.0;
   ++consdata->nvars;

   /* handle the new variable */
   SCIP_CALL( handleNewVariableSOS1(scip, cons, consdata, var, transformed) );

   return SCIP_OKAY;
}


/** deletes a variable of an SOS1 constraint */
static
SCIP_RETCODE deleteVarSOS1(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint */
   SCIP_CONSDATA*        consdata,           /**< constraint data */
   SCIP_EVENTHDLR*       eventhdlr,          /**< corresponding event handler */
   int                   pos                 /**< position of variable in array */
   )
{
   int j;

   assert( 0 <= pos && pos < consdata->nvars );

   /* remove lock of variable */
   SCIP_CALL( unlockVariableSOS1(scip, cons, consdata->vars[pos]) );

   /* drop events on variable */
   SCIP_CALL( SCIPdropVarEvent(scip, consdata->vars[pos], SCIP_EVENTTYPE_BOUNDCHANGED, eventhdlr, (SCIP_EVENTDATA*)consdata, -1) );

   /* delete variable - need to copy since order is important */
   for (j = pos; j < consdata->nvars-1; ++j)
   {
      consdata->vars[j] = consdata->vars[j+1]; /*lint !e679*/
      if ( consdata->weights != NULL )
         consdata->weights[j] = consdata->weights[j+1]; /*lint !e679*/
   }
   --consdata->nvars;

   return SCIP_OKAY;
}


/* ----------------------------- presolving --------------------------------------*/

/** extend a given clique in the conflict graph
 *
 *  Implementation of the Bron-Kerbosch Algorithm from the paper:
 *  Algorithm 457: Finding all Cliques of an Undirected Graph, Bron & Kerbosch, Commun. ACM, 1973
 */
static
SCIP_RETCODE extensionOperatorSOS1(
   SCIP*                 scip,               /**< SCIP pointer */
   SCIP_CONSHDLRDATA*    conshdlrdata,       /**< constraint handler data */
   SCIP_Bool**           adjacencymatrix,    /**< adjacencymatrix of the conflictgraph (only lower half filled) */
   SCIP_DIGRAPH*         vertexcliquegraph,  /**< graph that contains the information which cliques contain a given vertex
                                              *   vertices of variables = 0, ..., nsos1vars-1; vertices of cliques = nsos1vars, ..., nsos1vars+ncliques-1*/
   int                   nsos1vars,          /**< number of SOS1 variables */
   SCIP_CONS*            cons,               /**< constraint to be extended */
   SCIP_VAR**            vars,               /**< variables of extended clique */
   SCIP_Real*            weights,            /**< weights of extended clique */
   int**                 cliques,            /**< all cliques found so far */
   int*                  ncliques,           /**< number of clique found so far */
   int*                  cliquesizes,        /**< number of variables of current clique */
   int*                  newclique,          /**< clique we want to extended*/
   int*                  workingset,         /**< set of vertices that already served as extension and set of candidates that probably will lead to an extension */
   int                   nworkingset,        /**< length of array workingset */
   int                   nexts,              /**< number of vertices that already served as extension */
   int*                  maxextensions,      /**< maximal number of extensions */
   int*                  naddconss,          /**< number of added constraints */
   SCIP_Bool*            success             /**< pointer to store if at least one new clique was found */
   )
{
   int* workingsetnew;
   int nextsnew;
   int nworkingsetnew;
   int mincands;
   int btriter = 0; /* backtrack iterator */
   int selvertex;
   int selpos = -1;
   int fixvertex = -1;
   int i;
   int j;

   assert( scip != NULL );
   assert( conshdlrdata != NULL );
   assert( adjacencymatrix != NULL );
   assert( vertexcliquegraph != NULL );
   assert( cons != NULL );
   assert( cliques != NULL );
   assert( cliquesizes != NULL );
   assert( newclique != NULL );
   assert( workingset != NULL );
   assert( maxextensions != NULL );
   assert( naddconss != NULL );
   assert( success != NULL );

   *success = FALSE;

   mincands = nworkingset;

   SCIP_CALL( SCIPallocBufferArray(scip, &workingsetnew, nworkingset) );

#ifdef SCIP_DEBUG
   for (i = 0; i < nexts; ++i)
   {
      int vertex = workingset[i];
      for (j = nexts; j < nworkingset; ++j)
      {
         assert( isConnectedSOS1(adjacencymatrix, NULL, vertex, workingset[j]) );
      }
   }
#endif

   /* determine candidate with minimum number of disconnections */
   for (i = nexts; i < nworkingset && mincands != 0; ++i)
   {
      int pos = -1;
      int vertex = workingset[i];
      int cnt = 0;

      /* count disconnections */
      for (j = nexts; j < nworkingset && cnt < mincands; ++j)
      {
         if ( vertex != workingset[j] && ! isConnectedSOS1(adjacencymatrix, NULL, vertex, workingset[j]) )
         {
            cnt++;
            pos = j;
         }
      }

      /* check whether a new minimum was found */
      if ( cnt < mincands )
      {
         fixvertex = vertex;
         mincands = cnt;
         if ( i < nexts )
         {
            assert( pos >= 0 );
            selpos = pos;
         }
         else
         {
            selpos = i;
            btriter = 1;
         }
      }
   }

   /* If fixed point is initially chosen from candidates then number of disconnections will be preincreased by one. */

   /* Backtrackcycle */
   for (btriter = mincands + btriter; btriter >= 1; --btriter)
   {
      assert( selpos >= 0);
      assert( fixvertex >= 0);

      /* interchange */
      selvertex = workingset[selpos];
      workingset[selpos] = workingset[nexts];
      workingset[nexts] = selvertex;

      /* create new workingset */
      nextsnew = 0;
      for (j = 0 ; j < nexts; ++j)
      {
         if ( isConnectedSOS1(adjacencymatrix, NULL, selvertex, workingset[j]) )
            workingsetnew[nextsnew++] = workingset[j];
      }
      nworkingsetnew = nextsnew;
      for (j = nexts + 1; j < nworkingset; ++j)
      {
         if ( isConnectedSOS1(adjacencymatrix, NULL, selvertex, workingset[j]) )
            workingsetnew[nworkingsetnew++] = workingset[j];
      }

      newclique[cliquesizes[*ncliques]++] = selvertex;

      /* if we found a new clique */
      if ( nworkingsetnew == 0 )
      {
         char consname[SCIP_MAXSTRLEN];
         SCIP_CONSDATA* consdata;
         SCIP_CONS* newcons;
         int cliqueind;

         cliqueind = nsos1vars + *ncliques; /* index of clique in the vertex-clique graph */

         /* save new clique */
         SCIP_CALL( SCIPallocBufferArray(scip, &(cliques[*ncliques]), cliquesizes[*ncliques]) );
         for (j = 0 ; j < cliquesizes[*ncliques]; ++j)
         {
            vars[j] = nodeGetVarSOS1(conshdlrdata->conflictgraph, newclique[j]);
            weights[j] = j+1;
            cliques[*ncliques][j] = newclique[j];
         }

         SCIPsortInt(cliques[*ncliques], cliquesizes[*ncliques]);

         /* create new constraint */
         (void) SCIPsnprintf(consname, SCIP_MAXSTRLEN, "extsos1_%"SCIP_LONGINT_FORMAT, conshdlrdata->cntextsos1, conshdlrdata->cntextsos1);

         SCIP_CALL( SCIPcreateConsSOS1(scip, &newcons, consname, cliquesizes[*ncliques], vars, weights,
               SCIPconsIsInitial(cons), SCIPconsIsSeparated(cons),
               SCIPconsIsEnforced(cons), SCIPconsIsChecked(cons),
               SCIPconsIsPropagated(cons), SCIPconsIsLocal(cons),
               SCIPconsIsDynamic(cons),
               SCIPconsIsRemovable(cons), SCIPconsIsStickingAtNode(cons)) );

         consdata = SCIPconsGetData(newcons);

         /* add directed edges to the vertex-clique graph */
         for (j = 0; j < consdata->nvars; ++j)
            SCIP_CALL( SCIPdigraphAddArcSafe(vertexcliquegraph, cliques[*ncliques][j], cliqueind, NULL) );

         SCIP_CALL( SCIPaddCons(scip, newcons) );
         SCIP_CALL( SCIPreleaseCons(scip, &newcons) );

         ++(*naddconss);
         ++(conshdlrdata->cntextsos1);
         ++(*ncliques);
         cliquesizes[*ncliques] = cliquesizes[*ncliques-1]; /* cliquesizes[*ncliques] = size of newclique */

         *success = TRUE;

         --(*maxextensions);

         if ( *maxextensions <= 0 )
         {
            SCIPfreeBufferArray(scip, &workingsetnew);
            return SCIP_OKAY;
         }
      }
      else if ( nextsnew < nworkingsetnew ) /* else if the number of of candidates equals zero */
      {
         SCIP_CALL( extensionOperatorSOS1(scip, conshdlrdata, adjacencymatrix, vertexcliquegraph, nsos1vars, cons, vars, weights,
               cliques, ncliques, cliquesizes, newclique, workingsetnew, nworkingsetnew, nextsnew, maxextensions, naddconss, success) );

         if ( *maxextensions <= 0 )
         {
            SCIPfreeBufferArray(scip, &workingsetnew);
            return SCIP_OKAY;
         }
      }

      /* remove selvertex from clique */
      --cliquesizes[*ncliques];

      /* add selvertex to the set of vertices that already served as extension */
      ++nexts;

      if ( btriter > 1 )
      {
         /* select a candidate that is not connected to the fixed vertex */
         for (j = nexts; j < nworkingset; ++j)
         {
            if ( ! isConnectedSOS1(adjacencymatrix, NULL, fixvertex, workingset[j]) )
            {
               selpos = j;
               break;
            }
         }
      }
   }

   SCIPfreeBufferArray(scip, &workingsetnew);

   return SCIP_OKAY;
}


/** generates conflict graph that is induced by the variables of a linear constraint */
static
SCIP_RETCODE SCIPgenConflictgraphLinearCons(
   SCIP_CONSHDLR*        conshdlr,           /**< SOS1 constraint handler */
   SCIP_DIGRAPH*         conflictgraphlin,   /**< conflict graph of linear constraint (nodes: 1, ..., nlinvars) */
   SCIP_DIGRAPH*         conflictgraphorig,  /**< original conflict graph (nodes: 1, ..., nsos1vars) */
   SCIP_VAR**            linvars,            /**< linear variables in linear constraint */
   int                   nlinvars,           /**< number of linear variables in linear constraint */
   int*                  posinlinvars        /**< posinlinvars[i] = position (index) of SOS1 variable i in linear constraint,
                                              *   posinlinvars[i]= -1 if @p i is not a SOS1 variable or not a variable of the linear constraint */
   )
{
   int indexinsosvars;
   int indexinlinvars;
   int* succ;
   int nsucc;
   int v;
   int s;

   assert( conflictgraphlin != NULL );
   assert( conflictgraphorig != NULL );
   assert( linvars != NULL );
   assert( posinlinvars != NULL );

   for (v = 1; v < nlinvars; ++v) /* we start with v = 1, since "indexinlinvars < v" (see below) is never fulfilled for v = 0 */
   {
      indexinsosvars = varGetNodeSOS1(conshdlr, linvars[v]);

      /* if linvars[v] is contained in at least one SOS1 constraint */
      if ( indexinsosvars >= 0 )
      {
         succ = SCIPdigraphGetSuccessors(conflictgraphorig, indexinsosvars);
         nsucc = SCIPdigraphGetNSuccessors(conflictgraphorig, indexinsosvars);

         for (s = 0; s < nsucc; ++s)
         {
            assert( succ[s] >= 0 );
            indexinlinvars = posinlinvars[succ[s]];
            assert( indexinlinvars < nlinvars );

            if ( indexinlinvars >= 0 && indexinlinvars < v )
            {
               SCIP_CALL( SCIPdigraphAddArcSafe(conflictgraphlin, v, indexinlinvars, NULL) );
               SCIP_CALL( SCIPdigraphAddArcSafe(conflictgraphlin, indexinlinvars, v, NULL) );
            }
         }
      }
   }

   return SCIP_OKAY;
}


/** determine the common successors of the vertices from the considered clique */
static
SCIP_RETCODE CliqueGetCommonSuccessorsSOS1(
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler */
   SCIP_DIGRAPH*         conflictgraph,      /**< conflict graph */
   int*                  clique,             /**< current clique */
   SCIP_VAR**            vars,               /**< clique variables */
   int                   nvars,              /**< number of clique variables */
   int*                  comsucc,            /**< pointer to store common successors of clique vertices (size = nvars) */
   int*                  ncomsucc            /**< pointer to store number common successors of clique vertices */
   )
{
   int nsucc;
   int* succ;
   int ind;
   int k = 0;
   int v;
   int i;
   int j;

   assert( conflictgraph != NULL );
   assert( clique != NULL );
   assert( vars != NULL );
   assert( comsucc != NULL );
   assert( ncomsucc != NULL );

   *ncomsucc = 0;

   /* determine the common successors of the vertices from the considered clique */

   /* determine successors of variable var[0] that are not in the clique */
   ind =  varGetNodeSOS1(conshdlr, vars[0]);
   assert( ind != -1 );

   nsucc = SCIPdigraphGetNSuccessors(conflictgraph, ind);
   succ = SCIPdigraphGetSuccessors(conflictgraph, ind);

   for (j = 0; j < nvars; ++j)
   {
      for (i = k; i < nsucc; ++i)
      {
         if ( succ[i] > clique[j] )
         {
            k = i;
            break;
         }
         else if ( succ[i] == clique[j] )
         {
            k = i + 1;
            break;
         }
         else
            comsucc[(*ncomsucc)++] = succ[i];
      }
   }

   /* for all variables except the first one */
   for (v = 1; v < nvars; ++v)
   {
      int ncomsuccsave = 0;
      k = 0;

      ind =  varGetNodeSOS1(conshdlr, vars[v]);
      nsucc = SCIPdigraphGetNSuccessors(conflictgraph, ind);
      succ = SCIPdigraphGetSuccessors(conflictgraph, ind);

      /* determine successors that are in comsucc */
      for (j = 0; j < *ncomsucc; ++j)
      {
         for (i = k; i < nsucc; ++i)
         {
            if ( succ[i] > comsucc[j] )
            {
               k = i;
               break;
            }
            else if ( succ[i] == comsucc[j] )
            {
               comsucc[ncomsuccsave++] = succ[i];
               k = i + 1;
               break;
            }
         }
      }
      *ncomsucc = ncomsuccsave;
   }

   return SCIP_OKAY;
}


/* get nodes whose corresponding SOS1 variables are nonzero if an SOS1 variable of a given node is nonzero */
static
SCIP_RETCODE SCIPgetSOS1Implications(
   SCIP*                 scip,               /**< SCIP pointer */
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler */
   SCIP_VAR**            vars,               /**< problem and SOS1 variables */
   SCIP_DIGRAPH*         implgraph,          /**< implication graph (@j is successor of @p i if and only if \f$ x_i\not = 0 \Leftarrow x_j\not = 0\f$) */
   SCIP_HASHMAP*         implhash,           /**< hash map from variable to node in implication graph */
   SCIP_Bool*            implnodes,          /**< implnodes[i] = TRUE if the SOS1 variable corresponding to node i in the implication graph is implied to be nonzero */
   int                   node                /**< node of the implication graph */
   )
{
   int sos1node;

   assert( scip != NULL );
   assert( implgraph != NULL );
   assert( implnodes != NULL );
   assert( node >= 0 );
   assert( vars[node] != NULL );
   assert( (int) (size_t) SCIPhashmapGetImage(implhash, vars[node]) == node );

   /* get node of variable in the conflict graph (-1 if variable is no SOS1 variable) */
   sos1node = varGetNodeSOS1(conshdlr, vars[node]);
   if ( sos1node < 0 )
      return SCIP_OKAY;
   else
   {
      SCIP_SUCCDATA** succdatas;
      int* succ;
      int nsucc;
      int s;

      succdatas = (SCIP_SUCCDATA**) SCIPdigraphGetSuccessorsData(implgraph, node);
      nsucc = SCIPdigraphGetNSuccessors(implgraph, node);
      succ = SCIPdigraphGetSuccessors(implgraph, node);

      for (s = 0; s < nsucc; ++s)
      {
         SCIP_SUCCDATA* data;
         int succnode;
         succnode = succ[s];
         data = succdatas[s];
         sos1node = varGetNodeSOS1(conshdlr, vars[succnode]);

         /* if node is SOS1 and the corresponding variable is implied to be nonzero */
         assert( succdatas[s] != NULL );
         if ( sos1node >= 0 && ! implnodes[sos1node] && ( SCIPisFeasPositive(scip, data->lbimpl) || SCIPisFeasNegative(scip, data->ubimpl) ) )
         {
            assert( sos1node == succnode );
            implnodes[sos1node] = TRUE;
            SCIP_CALL( SCIPgetSOS1Implications(scip, conshdlr, vars, implgraph, implhash, implnodes, succnode) );
         }
      }
   }

   return SCIP_OKAY;
}


/* adds further edges to the conflict graph based on the information of an implication graph */
static
SCIP_RETCODE SCIPupdateConflictGraphSOS1(
   SCIP*                 scip,               /**< SCIP pointer */
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler */
   SCIP_DIGRAPH*         conflictgraph,      /**< conflict graph */
   SCIP_VAR**            totalvars,          /**< problem and SOS1 variables */
   SCIP_DIGRAPH*         implgraph,          /**< implication graph (@j is successor of @p i if and only if \f$ x_i\not = 0 \Leftarrow x_j\not = 0\f$) */
   SCIP_HASHMAP*         implhash,           /**< hash map from variable to node in implication graph */
   SCIP_Bool*            implnodes,          /**< implnodes[i] = TRUE if the SOS1 variable corresponding to a given node i in the implication graph is implied to be nonzero */
   SCIP_Bool**           adjacencymatrix,    /**< adjacencymatrix of the conflictgraph (only lower half filled) */
   int                   givennode,          /**< node of the conflict graph */
   int                   nonznode,           /**< node of the conflict graph that is implied to be nonzero if given node is nonzero */
   int*                  naddconss           /**< number of added SOS1 constraints */
   )
{
   SCIP_SUCCDATA** succdatas;
   int succnode;
   int* succ;
   int nsucc;
   int s;

   assert( nonznode >= 0 && nonznode < SCIPdigraphGetNNodes(conflictgraph) );

   nsucc = SCIPdigraphGetNSuccessors(conflictgraph, nonznode);
   succ = SCIPdigraphGetSuccessors(conflictgraph, nonznode);

   /* update conflict graph */
   for (s = 0; s < nsucc; ++s)
   {
      succnode = succ[s];

      if ( givennode != succnode )
      {
         if ( ! isConnectedSOS1(adjacencymatrix, NULL, givennode, succnode) )
         {
            char namesos[SCIP_MAXSTRLEN];
            SCIP_CONS* soscons = NULL;
            SCIP_VAR* var1;
            SCIP_VAR* var2;

            /* add arcs to the conflictgraph */
            SCIP_CALL( SCIPdigraphAddArc(conflictgraph, givennode, succnode, NULL) );
            SCIP_CALL( SCIPdigraphAddArc(conflictgraph, succnode, givennode, NULL) );

            /* update adjacencymatrix */
            if ( givennode > succnode )
               adjacencymatrix[givennode][succnode] = 1;
            else
               adjacencymatrix[succnode][givennode] = 1;

            var1 = nodeGetVarSOS1(conflictgraph, givennode);
            var2 = nodeGetVarSOS1(conflictgraph, succnode);

            /* create SOS1 constraint */
            (void) SCIPsnprintf(namesos, SCIP_MAXSTRLEN, "presolved_sos1_%s_%s", SCIPvarGetName(var1), SCIPvarGetName(var2) );
            SCIP_CALL( SCIPcreateConsSOS1(scip, &soscons, namesos, 0, NULL, NULL, TRUE, TRUE, TRUE, FALSE, TRUE,
                  TRUE, FALSE, FALSE, FALSE) );

            /* add variables to SOS1 constraint */
            SCIP_CALL( SCIPaddVarSOS1(scip, soscons, var1, 1.0) );
            SCIP_CALL( SCIPaddVarSOS1(scip, soscons, var2, 2.0) );

            /* add constraint */
            SCIP_CALL( SCIPaddCons(scip, soscons) );

            /* release constraint */
            SCIP_CALL( SCIPreleaseCons(scip, &soscons) );

            ++(*naddconss);
         }
      }
   }

   assert( nonznode == (int) (size_t) SCIPhashmapGetImage(implhash, nodeGetVarSOS1(conflictgraph, nonznode)) ); /* by constr.: nodes of SOS1 variables are equal for conflict graph and implication graph */
   succdatas = (SCIP_SUCCDATA**) SCIPdigraphGetSuccessorsData(implgraph, nonznode);
   nsucc = SCIPdigraphGetNSuccessors(implgraph, nonznode);
   succ = SCIPdigraphGetSuccessors(implgraph, nonznode);

   /* go further in implication graph */
   for (s = 0; s < nsucc; ++s)
   {
      SCIP_SUCCDATA* data;
      succnode = succ[s];
      data = succdatas[s];

      /* if node is SOS1 and not already known to be an implication node */
      if ( varGetNodeSOS1(conshdlr, totalvars[succnode]) >= 0 && ! implnodes[succnode] && ( SCIPisFeasPositive(scip, data->lbimpl) || SCIPisFeasNegative(scip, data->ubimpl) ) )
      {
         assert( succnode == (int) (size_t) SCIPhashmapGetImage(implhash, nodeGetVarSOS1(conflictgraph, succnode)) ); /* by constr.: nodes of SOS1 variables are equal for conflict graph and implication graph */
         implnodes[succnode] = TRUE; /* in order to avoid cycling */
         SCIP_CALL( SCIPupdateConflictGraphSOS1(scip, conshdlr, conflictgraph, totalvars, implgraph, implhash, implnodes, adjacencymatrix, givennode, succnode, naddconss) );
      }
   }

   return SCIP_OKAY;
}

/** perform one presolving round for a single SOS1 constraint
 *
 *  We perform the following presolving steps.
 *
 *  - If the bounds of some variable force it to be nonzero, we can
 *    fix all other variables to zero and remove the SOS1 constraints
 *    that contain it.
 *  - If a variable is fixed to zero, we can remove the variable.
 *  - If a variable appears twice, it can be fixed to 0.
 *  - We substitute appregated variables.
 */
static
SCIP_RETCODE presolRoundConsSOS1(
   SCIP*                 scip,               /**< SCIP pointer */
   SCIP_CONS*            cons,               /**< constraint */
   SCIP_CONSDATA*        consdata,           /**< constraint data */
   SCIP_EVENTHDLR*       eventhdlr,          /**< event handler */
   SCIP_Bool*            cutoff,             /**< whether a cutoff happened */
   SCIP_Bool*            success,            /**< whether we performed a successful reduction */
   int*                  ndelconss,          /**< number of deleted constraints */
   int*                  nupgdconss,         /**< number of upgraded constraints */
   int*                  nfixedvars,         /**< number of fixed variables */
   int*                  nremovedvars        /**< number of variables removed */
   )
{
   SCIP_VAR** vars;
   SCIP_Bool allvarsbinary;
   SCIP_Bool infeasible;
   SCIP_Bool fixed;
   int nfixednonzeros;
   int lastFixedNonzero;
   int j;

   assert( scip != NULL );
   assert( cons != NULL );
   assert( consdata != NULL );
   assert( eventhdlr != NULL );
   assert( cutoff != NULL );
   assert( success != NULL );
   assert( ndelconss != NULL );
   assert( nfixedvars != NULL );
   assert( nremovedvars != NULL );

   *cutoff = FALSE;
   *success = FALSE;

   SCIPdebugMessage("Presolving SOS1 constraint <%s>.\n", SCIPconsGetName(cons) );

   j = 0;
   nfixednonzeros = 0;
   lastFixedNonzero = -1;
   allvarsbinary = TRUE;
   vars = consdata->vars;

   /* check for variables fixed to 0 and bounds that fix a variable to be nonzero */
   while ( j < consdata->nvars )
   {
      int l;
      SCIP_VAR* var;
      SCIP_Real lb;
      SCIP_Real ub;
      SCIP_Real scalar;
      SCIP_Real constant;

      scalar = 1.0;
      constant = 0.0;

      /* check for aggregation: if the constant is zero the variable is zero iff the aggregated
       * variable is 0 */
      var = vars[j];
      SCIP_CALL( SCIPgetProbvarSum(scip, &var, &scalar, &constant) );

      /* if constant is zero and we get a different variable, substitute variable */
      if ( SCIPisZero(scip, constant) && ! SCIPisZero(scip, scalar) && var != vars[j] )
      {
         SCIPdebugMessage("substituted variable <%s> by <%s>.\n", SCIPvarGetName(vars[j]), SCIPvarGetName(var));
         SCIP_CALL( SCIPdropVarEvent(scip, consdata->vars[j], SCIP_EVENTTYPE_BOUNDCHANGED, eventhdlr, (SCIP_EVENTDATA*)consdata, -1) );
         SCIP_CALL( SCIPcatchVarEvent(scip, var, SCIP_EVENTTYPE_BOUNDCHANGED, eventhdlr, (SCIP_EVENTDATA*)consdata, NULL) );

         /* change the rounding locks */
         SCIP_CALL( unlockVariableSOS1(scip, cons, consdata->vars[j]) );
         SCIP_CALL( lockVariableSOS1(scip, cons, var) );

         vars[j] = var;
      }

      /* check whether the variable appears again later */
      for (l = j+1; l < consdata->nvars; ++l)
      {
         /* if variable appeared before, we can fix it to 0 and remove it */
         if ( vars[j] == vars[l] )
         {
            SCIPdebugMessage("variable <%s> appears twice in constraint, fixing it to 0.\n", SCIPvarGetName(vars[j]));
            SCIP_CALL( SCIPfixVar(scip, vars[j], 0.0, &infeasible, &fixed) );

            if ( infeasible )
            {
               *cutoff = TRUE;
               return SCIP_OKAY;
            }
            if ( fixed )
               ++(*nfixedvars);
         }
      }

      /* get bounds */
      lb = SCIPvarGetLbLocal(vars[j]);
      ub = SCIPvarGetUbLocal(vars[j]);

      /* if the variable if fixed to nonzero */
      if ( SCIPisFeasPositive(scip, lb) || SCIPisFeasNegative(scip, ub) )
      {
         ++nfixednonzeros;
         lastFixedNonzero = j;
      }

      /* if the variable is fixed to 0 */
      if ( SCIPisFeasZero(scip, lb) && SCIPisFeasZero(scip, ub) )
      {
         SCIPdebugMessage("deleting variable <%s> fixed to 0.\n", SCIPvarGetName(vars[j]));
         SCIP_CALL( deleteVarSOS1(scip, cons, consdata, eventhdlr, j) );
         ++(*nremovedvars);
      }
      else
      {
         /* check whether all variables are binary */
         if ( ! SCIPvarIsBinary(vars[j]) )
            allvarsbinary = FALSE;

         ++j;
      }
   }

   /* if the number of variables is less than 2 */
   if ( consdata->nvars < 2 )
   {
      SCIPdebugMessage("Deleting SOS1 constraint <%s> with < 2 variables.\n", SCIPconsGetName(cons));

      /* delete constraint */
      assert( ! SCIPconsIsModifiable(cons) );
      SCIP_CALL( SCIPdelCons(scip, cons) );
      ++(*ndelconss);
      *success = TRUE;
      return SCIP_OKAY;
   }

   /* if more than one variable are fixed to be nonzero, we are infeasible */
   if ( nfixednonzeros > 1 )
   {
      SCIPdebugMessage("The problem is infeasible: more than one variable has bounds that keep it from being 0.\n");
      assert( lastFixedNonzero >= 0 );
      *cutoff = TRUE;
      return SCIP_OKAY;
   }

   /* if there is exactly one fixed nonzero variable */
   if ( nfixednonzeros == 1 )
   {
      assert( lastFixedNonzero >= 0 );

      /* fix all other variables to zero */
      for (j = 0; j < consdata->nvars; ++j)
      {
         if ( j != lastFixedNonzero )
         {
            SCIP_CALL( SCIPfixVar(scip, vars[j], 0.0, &infeasible, &fixed) );
            assert( ! infeasible );
            if ( fixed )
               ++(*nfixedvars);
         }
      }

      SCIPdebugMessage("Deleting redundant SOS1 constraint <%s> with one variable.\n", SCIPconsGetName(cons));

      /* delete original constraint */
      assert( ! SCIPconsIsModifiable(cons) );
      SCIP_CALL( SCIPdelCons(scip, cons) );
      ++(*ndelconss);
      *success = TRUE;
   }
   /* note: there is no need to update consdata->nfixednonzeros, since the constraint is deleted as soon nfixednonzeros > 0. */
   else
   {
      /* if all variables are binary create a set packing constraint */
      if ( allvarsbinary )
      {
         SCIP_CONS* setpackcons;

         /* create, add, and release the logicor constraint */
         SCIP_CALL( SCIPcreateConsSetpack(scip, &setpackcons, SCIPconsGetName(cons), consdata->nvars, consdata->vars,
               SCIPconsIsInitial(cons), SCIPconsIsSeparated(cons), SCIPconsIsEnforced(cons), SCIPconsIsChecked(cons),
               SCIPconsIsPropagated(cons), SCIPconsIsLocal(cons), SCIPconsIsModifiable(cons), SCIPconsIsDynamic(cons), 
               SCIPconsIsRemovable(cons), SCIPconsIsStickingAtNode(cons)) );
         SCIP_CALL( SCIPaddCons(scip, setpackcons) );
         SCIP_CALL( SCIPreleaseCons(scip, &setpackcons) );

         SCIPdebugMessage("Upgrading SOS1 constraint <%s> to set packing constraint.\n", SCIPconsGetName(cons));

         /* remove the SOS1 constraint globally */
         assert( ! SCIPconsIsModifiable(cons) );
         SCIP_CALL( SCIPdelCons(scip, cons) );
         ++(*nupgdconss);
         *success = TRUE;
      }
   }

   return SCIP_OKAY;
}



/** perform one presolving round for all SOS1 constraints
 *
 *  We perform the following presolving steps.
 *
 *  - If the bounds of some variable force it to be nonzero, we can
 *    fix all other variables to zero and remove the SOS1 constraints
 *    that contain it.
 *  - If a variable is fixed to zero, we can remove the variable.
 *  - If a variable appears twice, it can be fixed to 0.
 *  - We substitute appregated variables.
 *  - Search for larger SOS1 constraints in the conflict graph
 *  - Remove redundant SOS1 constraints
 */
static
SCIP_RETCODE presolRoundConssSOS1(
   SCIP*                 scip,               /**< SCIP pointer */
   SCIP_EVENTHDLR*       eventhdlr,          /**< event handler */
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler */
   SCIP_CONSHDLRDATA*    conshdlrdata,       /**< constraint handler data */
   SCIP_DIGRAPH*         conflictgraph,      /**< conflict graph */
   SCIP_Bool**           adjacencymatrix,    /**< adjacency matrix of conflict graph */
   SCIP_CONS**           conss,              /**< SOS1 constraints */
   int                   nconss,             /**< number of SOS1 constraints */
   int                   nsos1vars,          /**< number of SOS1 variables */
   int*                  naddconss,          /**< number of added constraints */
   int*                  ndelconss,          /**< number of deleted constraints */
   int*                  nupgdconss,         /**< number of upgraded constraints */
   int*                  nfixedvars,         /**< number of fixed variables */
   int*                  nremovedvars,       /**< number of variables removed */
   SCIP_RESULT*          result              /**< result */
   )
{
   SCIP_DIGRAPH* vertexcliquegraph;
   SCIP_VAR** consvars;
   SCIP_Real* consweights;
   int** cliques = NULL;
   int ncliques = 0;
   int* cliquesizes = NULL;
   int* newclique = NULL;
   int* indconss = NULL;
   int* lengthconss = NULL;
   int* comsucc = NULL;
   int csize;
   int iter;
   int c;
   int i;

   /* create digraph whose nodes represent variables and cliques in the conflict graph */
   csize = MAX(1, conshdlrdata->maxextensions) * nconss;
   SCIP_CALL( SCIPdigraphCreate(&vertexcliquegraph, nsos1vars + csize) );

   /* allocate buffer arrays */
   SCIP_CALL( SCIPallocBufferArray(scip, &consvars, nsos1vars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &consweights, nsos1vars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &cliquesizes, csize) );
   SCIP_CALL( SCIPallocBufferArray(scip, &newclique, nsos1vars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &indconss, csize) );
   SCIP_CALL( SCIPallocBufferArray(scip, &lengthconss, csize) );
   SCIP_CALL( SCIPallocBufferArray(scip, &comsucc, MAX(nsos1vars, csize)) );
   SCIP_CALL( SCIPallocBufferArray(scip, &cliques, csize) );

   /* get constraint indices and sort them in descending order of their lengths */
   for (c = 0; c < nconss; ++c)
   {
      SCIP_CONSDATA* consdata;

      consdata = SCIPconsGetData(conss[c]);
      assert( consdata != NULL );

      indconss[c] = c;
      lengthconss[c] = consdata->nvars;
   }
   SCIPsortDownIntInt(lengthconss, indconss, nconss);

   /* check each constraint */
   for (iter = 0; iter < nconss; ++iter)
   {
      SCIP_CONSDATA* consdata;
      SCIP_CONS* cons;
      SCIP_Bool success;
      SCIP_Bool cutoff;
      int savendelconss;

      SCIP_VAR** vars;
      int nvars;

      c = indconss[iter];

      assert( conss != NULL );
      assert( conss[c] != NULL );
      cons = conss[c];
      consdata = SCIPconsGetData(cons);

      assert( consdata != NULL );
      assert( consdata->nvars >= 0 );
      assert( consdata->nvars <= consdata->maxvars );
      assert( ! SCIPconsIsModifiable(cons) );

      savendelconss = *ndelconss;

      /* perform one presolving round for SOS1 constraint */
      SCIP_CALL( presolRoundConsSOS1(scip, cons, consdata, eventhdlr, &cutoff, &success, ndelconss, nupgdconss, nfixedvars, nremovedvars) );

      if ( cutoff )
      {
         *result = SCIP_CUTOFF;
         break;
      }

      if ( *ndelconss > savendelconss )
      {
         *result = SCIP_SUCCESS;
         continue;
      }

      if ( success )
         *result = SCIP_SUCCESS;

      /* get number of variables of constraint */
      nvars = consdata->nvars;

      /* get variables of constraint */
      vars = consdata->vars;

      if ( nvars > 1 && conshdlrdata->maxextensions != 0 && ( nconss <= conshdlrdata->maxconsdelayext || conshdlrdata->maxconsdelayext == -1 ) )
      {
         SCIP_Bool extended = FALSE;
         int cliquesize = 0;
         int ncomsucc = 0;
         int varprobind;
         int k;
         int j;

         /* get clique and size of clique */
         for (j = 0; j < nvars; ++j)
         {
            varprobind = varGetNodeSOS1(conshdlr, vars[j]);

            if ( varprobind >= 0 )
               newclique[cliquesize++] = varprobind;
         }

         if ( cliquesize > 1 )
         {
            cliquesizes[ncliques] = cliquesize;

            /* sort clique vertices */
            SCIPsortInt(newclique, cliquesizes[ncliques]);

            /* check if clique is contained in an already known clique */
            if ( ncliques > 0 )
            {
               int* succ;
               int nsucc;
               int v;

               varprobind = newclique[0];
               ncomsucc = SCIPdigraphGetNSuccessors(vertexcliquegraph, varprobind);
               succ = SCIPdigraphGetSuccessors(vertexcliquegraph, varprobind);

               for (j = 0; j < ncomsucc; ++j)
               {
                  /* successors should have been sorted in a former step of the algorithm */
                  assert( j == 0 || succ[j] > succ[j-1] );
                  comsucc[j] = succ[j];
               }

               for (v = 1; v < cliquesize; ++v)
               {
                  int ncomsuccsave = 0;
                  k = 0;

                  varprobind = newclique[v];
                  nsucc = SCIPdigraphGetNSuccessors(vertexcliquegraph, varprobind);
                  succ = SCIPdigraphGetSuccessors(vertexcliquegraph, varprobind);

                  /* determine successors that are in comsucc */
                  for (j = 0; j < ncomsucc; ++j)
                  {
                     for (i = k; i < nsucc; ++i)
                     {
                        assert( i == 0 || succ[i] > succ[i-1] );

                        if ( succ[i] > comsucc[j] )
                        {
                           k = i;
                           break;
                        }
                        else if ( succ[i] == comsucc[j] )
                        {
                           comsucc[ncomsuccsave++] = succ[i];
                           k = i + 1;
                           break;
                        }
                     }
                  }
                  ncomsucc = ncomsuccsave;
               }
            }

            /* if constraint is redundand then delete it */
            if ( ncomsucc > 0 )
            {
               assert( ! SCIPconsIsModifiable(cons) );
               SCIP_CALL( SCIPdelCons(scip, cons) );
               ++(*ndelconss);
               *result = SCIP_SUCCESS;
               continue;
            }

            if ( conshdlrdata->maxextensions != 0 )
            {
               int maxextensions;
               ncomsucc = 0;

               /* determine the common successors of the vertices from the considered clique */
               SCIP_CALL( CliqueGetCommonSuccessorsSOS1(conshdlr, conflictgraph, newclique, vars, nvars, comsucc, &ncomsucc) );

               /* find extensions for the clique */
               maxextensions = conshdlrdata->maxextensions;

               SCIP_CALL( extensionOperatorSOS1(scip, conshdlrdata, adjacencymatrix, vertexcliquegraph, nsos1vars, cons, consvars, consweights,
                     cliques, &ncliques, cliquesizes, newclique, comsucc, ncomsucc, 0, &maxextensions, naddconss, &extended) );
            }

            /* if an extension was found for the current clique then free the old SOS1 constraint */
            if ( extended )
            {
               assert( ! SCIPconsIsModifiable(cons) );
               SCIP_CALL( SCIPdelCons(scip, cons) );
               ++(*ndelconss);
               *result = SCIP_SUCCESS;
            }
            else /* if we keep the constraint */
            {
               int cliqueind;

               cliqueind = nsos1vars + ncliques; /* index of clique in vertex-clique graph */

               /* add directed edges to the vertex-clique graph */
               SCIP_CALL( SCIPallocBufferArray(scip, &cliques[ncliques], cliquesize) );
               for (j = 0; j < cliquesize; ++j)
               {
                  cliques[ncliques][j] = newclique[j];
                  SCIP_CALL( SCIPdigraphAddArcSafe(vertexcliquegraph, cliques[ncliques][j], cliqueind, NULL) );
               }

               /* update number of maximal cliques */
               ++ncliques;
            }
         }
      }
   }

   /* free buffer arrays */
   for (c = ncliques-1; c >= 0; --c)
      SCIPfreeBufferArrayNull(scip, &cliques[c]);
   SCIPfreeBufferArrayNull(scip, &cliques);
   SCIPfreeBufferArrayNull(scip, &comsucc);
   SCIPfreeBufferArrayNull(scip, &lengthconss);
   SCIPfreeBufferArrayNull(scip, &indconss);
   SCIPfreeBufferArrayNull(scip, &newclique);
   SCIPfreeBufferArrayNull(scip, &cliquesizes);
   SCIPfreeBufferArrayNull(scip, &consweights);
   SCIPfreeBufferArrayNull(scip, &consvars);
   SCIPdigraphFree(&vertexcliquegraph);

   return SCIP_OKAY;
}


/* returns whether node is implied to be zero; this information is taken from the input array 'implnodes' */
static
SCIP_Bool isImpliedZero(
   SCIP_DIGRAPH*         conflictgraph,      /**< conflict graph */
   SCIP_Bool*            implnodes,          /**< implnodes[i] = TRUE if the SOS1 variable corresponding to node i in the implication graph is implied to be nonzero */
   int                   node                /**< node of the conflict graph (or -1) */
   )
{
   int* succ;
   int nsucc;
   int s;

   if ( node >= 0 )
   {
      nsucc = SCIPdigraphGetNSuccessors(conflictgraph, node);
      succ = SCIPdigraphGetSuccessors(conflictgraph, node);

      /* check whether any successor is implied to be nonzero */
      for (s = 0; s < nsucc; ++s)
      {
         if ( implnodes[succ[s]] )
            return TRUE;
      }
   }

   return FALSE;
}


/* updates arc data of implication graph */
static
SCIP_Bool updateArcData(
   SCIP*                 scip,               /**< SCIP pointer */
   SCIP_DIGRAPH*         implgraph,          /**< implication graph */
   SCIP_HASHMAP*         implhash,           /**< hash map from variable to node in implication graph */
   SCIP_VAR**            totalvars,          /**< problem and SOS1 variables */
   SCIP_VAR*             varv,               /**< variable that is assumed to be nonzero */
   SCIP_VAR*             varw,               /**< implication variable */
   SCIP_Real             lb,                 /**< old lower bound of \f$x_w\f$ */
   SCIP_Real             ub,                 /**< old upper bound of \f$x_w\f$ */
   SCIP_Real             newbound,           /**< new bound of \f$x_w\f$ */
   SCIP_Bool             lower,              /**< whether to consider lower bound implication (otherwise upper bound) */
   int*                  nchgbds,            /**< pointer to store number of changed bounds */
   SCIP_Bool*            update              /**< pointer to store whether implication graph has been updated */
   )
{
   SCIP_SUCCDATA** succdatas;
   SCIP_SUCCDATA* data = NULL;
   int nsucc;
   int* succ;
   int indv;
   int indw;
   int s;

   assert( scip != NULL );
   assert( implgraph != NULL );
   assert( implhash != NULL );
   assert( totalvars != NULL );
   assert( varv != NULL );
   assert( varw != NULL );

   /* if x_v != 0 turns out to be infeasible then fix x_v = 0 */
   if ( ( lower && SCIPisFeasLT(scip, ub, newbound) ) || ( ! lower && SCIPisFeasGT(scip, lb, newbound) ) )
   {
      SCIP_Bool infeasible;
      SCIP_Bool tightened1;
      SCIP_Bool tightened2;

      SCIP_CALL( SCIPtightenVarLb(scip, varv, 0.0, FALSE, &infeasible, &tightened1) );
      assert( !infeasible );
      SCIP_CALL( SCIPtightenVarUb(scip, varv, 0.0, FALSE, &infeasible, &tightened2) );
      assert( !infeasible );

      if ( tightened1 || tightened2 )
      {
         SCIPdebugMessage("fixed variable %s from lb = %f and ub = %f to 0.0 \n", SCIPvarGetName(varv), lb, ub);
         ++(*nchgbds);
      }
   }

   /* get successor information */
   indv = (int) (size_t) SCIPhashmapGetImage(implhash, varv); /* get index of x_v in implication graph */
   assert( (int) (size_t) SCIPhashmapGetImage(implhash, totalvars[indv]) == indv );
   succdatas = (SCIP_SUCCDATA**) SCIPdigraphGetSuccessorsData(implgraph, indv);
   nsucc = SCIPdigraphGetNSuccessors(implgraph, indv);
   succ = SCIPdigraphGetSuccessors(implgraph, indv);

   /* search for nodew in existing successors. If this is the case then check whether the lower implication bound may be updated ... */
   indw = (int) (size_t) SCIPhashmapGetImage(implhash, varw);
   assert( (int) (size_t) SCIPhashmapGetImage(implhash, totalvars[indw]) == indw );
   for (s = 0; s < nsucc; ++s)
   {
      if ( succ[s] == indw )
      {
         data = succdatas[s];
         assert( data != NULL );
         if ( lower && SCIPisFeasLT(scip, data->lbimpl, newbound) )
         {
            data->lbimpl = newbound;
            *update = TRUE;
            SCIPdebugMessage("updated to implication %s != 0 -> %s >= %f\n", SCIPvarGetName(varv), SCIPvarGetName(varw), newbound);
         }
         else if ( ! lower && SCIPisFeasGT(scip, data->ubimpl, newbound) )
         {
            data->ubimpl = newbound;
            *update = TRUE;
            SCIPdebugMessage("updated to implication %s != 0 -> %s >= %f\n", SCIPvarGetName(varv), SCIPvarGetName(varw), newbound);
         }
         break;
      }
   }

   /* ..., otherwise if there does not exist an arc between indv and indw already, then create one and add implication */
   if ( s == nsucc )
   {
      assert( data == NULL );
      SCIP_CALL( SCIPallocMemory(scip, &data) );
      if ( lower )
      {
         data->lbimpl = newbound;
         data->ubimpl = ub;
         SCIPdebugMessage("add implication %s != 0 -> %s >= %f\n", SCIPvarGetName(varv), SCIPvarGetName(varw), newbound);
      }
      else
      {
         data->lbimpl = lb;
         data->ubimpl = newbound;
         SCIPdebugMessage("add implication %s != 0 -> %s <= %f\n", SCIPvarGetName(varv), SCIPvarGetName(varw), newbound);
      }
      SCIP_CALL( SCIPdigraphAddArc(implgraph, indv, indw, (void*) data) );
      *update = TRUE;
   }

   return SCIP_OKAY;
}


/* assume some variable is nonzero. If this implies that some other variable is also nonzero, then store this information in an implication graph */
static
SCIP_RETCODE updateImplicationGraphSOS1(
   SCIP*                 scip,               /**< SCIP pointer */
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler */
   SCIP_DIGRAPH*         conflictgraph,      /**< conflict graph */
   SCIP_Bool**           adjacencymatrix,    /**< adjacency matrix of conflict graph (lower half) */
   SCIP_DIGRAPH*         implgraph,          /**< implication graph (@p j is successor of @p i if and only if \f$ x_i\not = 0 \Leftarrow x_j\not = 0\f$) */
   SCIP_HASHMAP*         implhash,           /**< hash map from variable to node in implication graph */
   SCIP_Bool*            implnodes,          /**< implnodes[i] = TRUE if the SOS1 variable corresponding to node i in the implication graph is implied to be nonzero */
   SCIP_VAR**            totalvars,          /**< problem and SOS1 variables */
   int**                 cliquecovers,       /**< clique covers of linear constraint */
   int*                  cliquecoversizes,   /**< size of clique covers */
   int*                  varincover,         /**< array with varincover[i] = cover of SOS1 index @p i */
   SCIP_VAR**            vars,               /**< variables to be checked */
   SCIP_Real*            coefs,              /**< coefficients of variables in linear constraint */
   int                   nvars,              /**< number of variables to be checked */
   SCIP_Real*            bounds,             /**< bounds of variables */
   SCIP_VAR*             var,                /**< variable that is assumed to be nonzero */
   SCIP_Real             bound,              /**< bound of variable */
   SCIP_Real             boundnonzero,       /**< bound of variable if it is known to be nonzero if infinity values are not summarized */
   int                   ninftynonzero,      /**< number of times infinity/-infinity has to be summarized to boundnonzero */
   SCIP_Bool             lower,              /**< TRUE if lower bounds are consideres; FALSE for upper bounds */
   int*                  nchgbds,            /**< pointer to store number of changed bounds */
   SCIP_Bool*            update              /**< pointer to store whether implication graph has been updated */
   )
{
   int nodev;
   int w;

   *update = FALSE;
   nodev = varGetNodeSOS1(conshdlr, var); /* possibly -1 if var is not involved in an SOS1 constraint */

   /* update implication graph if possible */

   /* if nodev is an index of an SOS1 variable and at least one lower bound of a variable that is not x_v is infinity */
   if ( nodev >= 0 && ! SCIPisInfinity(scip, REALABS(bound) ) && ninftynonzero <= 1 )
   {
      /* for every variable x_w: compute upper bound of a_w * x_w if x_v is known to be nonzero */
      for (w = 0; w < nvars; ++w)
      {
         int newninftynonzero;
         SCIP_Bool implinfty = FALSE;
         int nodew;

         /* get node of x_w in conflict graph: nodew = -1 if it is no SOS1 variable */
         nodew = varGetNodeSOS1(conshdlr, vars[w]);

         newninftynonzero = ninftynonzero;

         /* variable should not be fixed to be already zero (note x_v is fixed to be nonzero by assumption) */
         if ( nodew < 0 || ( nodev != nodew && ! isConnectedSOS1(adjacencymatrix, NULL, nodev, nodew) && ! isImpliedZero(conflictgraph, implnodes, nodew) ) )
         {
            SCIP_Real implbound;
            SCIP_Bool implcoverw;
            int nodecliq;
            int indcliq;
            int ind;
            int j;

            /* boundnonzero is the bound of x_v if x_v is nonzero
               we use this information to get a bound of x_w if x_v is nonzero;
               therefore, we have to perform some recomputations */
            implbound = boundnonzero - bound;
            ind = varincover[w];
            assert( cliquecoversizes[ind] > 0 );

            implcoverw = FALSE;
            for (j = 0; j < cliquecoversizes[ind]; ++j)
            {
               indcliq = cliquecovers[ind][j];
               assert( 0 <= indcliq && indcliq < nvars );

               nodecliq = varGetNodeSOS1(conshdlr, vars[indcliq]); /* possibly -1 if variable is not involved in an SOS1 constraint */

               /* if nodecliq is not a member of an SOS1 constraint or the variable corresponding to nodecliq is not implied to be zero if x_v != 0  */
               if ( nodecliq < 0 || (! isConnectedSOS1(adjacencymatrix, NULL, nodev, nodecliq) && ! isImpliedZero(conflictgraph, implnodes, nodecliq) ) )
               {
                  if ( indcliq == w )
                  {
                     if ( ! SCIPisInfinity(scip, REALABS(bounds[w])) )
                        implbound += bounds[w];
                     else
                        --newninftynonzero;
                     implcoverw = TRUE;
                  }
                  else if ( implcoverw )
                  {
                     if ( SCIPisInfinity(scip, REALABS( bounds[indcliq] )) )
                        implinfty = TRUE;
                     else
                        implbound -= bounds[indcliq];
                     break;
                  }
                  else
                  {
                     if ( SCIPisInfinity(scip, REALABS( bounds[indcliq] ) ) )
                        implinfty = TRUE;
                     break;
                  }
               }
            }

            /* check whether x_v != 0 implies a bound change of x_w */
            if ( ! implinfty && newninftynonzero == 0 )
            {
               SCIP_Real newbound;
               SCIP_Real coef;
               SCIP_Real lb;
               SCIP_Real ub;

               lb = SCIPvarGetLbLocal(vars[w]);
               ub = SCIPvarGetUbLocal(vars[w]);
               coef = coefs[w];
               assert( ! SCIPisFeasZero(scip, coef) );
               newbound = implbound / coef;

               /* check if an implication can be added/updated or assumption x_v != 0 is infeasible */
               if ( lower )
               {
                  if ( SCIPisFeasPositive(scip, coef) && SCIPisFeasLT(scip, lb, newbound) )
                  {
                     SCIP_CALL( updateArcData(scip, implgraph, implhash, totalvars, var, vars[w], lb, ub, newbound, TRUE, nchgbds, update) );
                  }
                  else if ( SCIPisFeasNegative(scip, coef) && SCIPisFeasGT(scip, ub, newbound) )
                  {
                     SCIP_CALL( updateArcData(scip, implgraph, implhash, totalvars, var, vars[w], lb, ub, newbound, FALSE, nchgbds, update) );
                  }
               }
               else
               {
                  if ( SCIPisFeasPositive(scip, coef) && SCIPisFeasGT(scip, ub, newbound) )
                  {
                     SCIP_CALL( updateArcData(scip, implgraph, implhash, totalvars, var, vars[w], lb, ub, newbound, FALSE, nchgbds, update) );
                  }
                  else if ( SCIPisFeasNegative(scip, coef) && SCIPisFeasLT(scip, lb, newbound) )
                  {
                     SCIP_CALL( updateArcData(scip, implgraph, implhash, totalvars, var, vars[w], lb, ub, newbound, TRUE, nchgbds, update) );
                  }
               }
            }
         }
      }
   }

   return SCIP_OKAY;
}


/** for a given index @p v search for a clique in the conflict graph restricted to the variables of a linear constraint that covers @p v
 *  and has an an empty intersection with already computed clique covers
 */
static
SCIP_RETCODE SCIPcomputeVarsCoverSOS1(
   SCIP*                 scip,               /**< SCIP pointer */
   SCIP_DIGRAPH*         conflictgraphroot,  /**< conflict graph of the root node (nodes: 1, ..., @p nsos1vars) */
   SCIP_DIGRAPH*         conflictgraphlin,   /**< conflict graph of linear constraint (nodes: 1, ..., @p nlinvars) */
   SCIP_VAR**            linvars,            /**< variables in linear constraint */
   SCIP_Bool*            coveredvars,        /**< states which variables of the linear constraint are currently covered by a clique */
   int*                  clique,             /**< array to store new clique in cover */
   int*                  cliquesize,         /**< pointer to store the size of @p clique */
   int                   v,                  /**< position of variable in linear constraint that should be covered */
   SCIP_Bool             considersolvals     /**< TRUE if largest auxiliary bigM values of variables should be prefered */
   )
{
   int nsucc;
   int s;

   assert( conflictgraphlin != NULL );
   assert( linvars != NULL );
   assert( coveredvars != NULL );
   assert( clique != NULL );
   assert( cliquesize != NULL );

   assert( ! coveredvars[v] );  /* we should produce a new clique */

   /* add index 'v' to the clique cover */
   clique[0] = v;
   *cliquesize = 1;

   nsucc = SCIPdigraphGetNSuccessors(conflictgraphlin, v);
   if ( nsucc > 0 )
   {
      int* extensions;
      int nextensions = 0;
      int nextensionsnew;
      int succnode;
      int* succ;

      /* allocate buffer array */
      SCIP_CALL( SCIPallocBufferArray(scip, &extensions, nsucc) );

      succ = SCIPdigraphGetSuccessors(conflictgraphlin, v);

      /* compute possible extensions for the clique cover */
      for (s = 0; s < nsucc; ++s)
      {
         succnode = succ[s];
         if ( ! coveredvars[succnode] )
            extensions[nextensions++] = succ[s];
      }

      /* while there exist possible extensions for the clique cover */
      while ( nextensions > 0 )
      {
         int bestindex = -1;

         if ( considersolvals )
         {
            SCIP_Real bestbigMval;
            SCIP_Real bigMval;

            bestbigMval = -SCIPinfinity(scip);

            /* search for the extension with the largest absolute value of its LP relaxation solution value */
            for (s = 0; s < nextensions; ++s)
            {
               bigMval = SCIPnodeGetSolvalBinaryBigMSOS1(scip, conflictgraphroot, NULL, extensions[s]);
               if ( SCIPisFeasLT(scip, bestbigMval, bigMval) )
               {
                  bestbigMval = bigMval;
                  bestindex = extensions[s];
               }
            }
         }
         else
            bestindex = extensions[0];

         assert( bestindex != -1 );

         /* add bestindex to the clique cover */
         clique[(*cliquesize)++] = bestindex;

         /* compute new 'extensions' array */
         nextensionsnew = 0;
         for (s = 0; s < nextensions; ++s)
         {
            if ( s != bestindex && isConnectedSOS1(NULL, conflictgraphlin, bestindex, extensions[s]) )
               extensions[nextensionsnew++] = extensions[s];
         }
         nextensions = nextensionsnew;
      }

      /* free buffer array */
      SCIPfreeBufferArray(scip, &extensions);
   }

   /* mark covered indices */
   for (s = 0; s < *cliquesize; ++s)
   {
      int ind = clique[s];

      assert( ! coveredvars[ind] );
      coveredvars[ind] = TRUE;
   }

   return SCIP_OKAY;
}


/* try to tighten upper and lower bounds for variables */
static
SCIP_RETCODE tightenVarsBoundsSOS1(
   SCIP*                 scip,               /**< SCIP pointer */
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler */
   SCIP_DIGRAPH*         conflictgraph,      /**< conflict graph */
   SCIP_DIGRAPH*         implgraph,          /**< implication graph (@j is successor of @p i if and only if \f$ x_i\not = 0 \f$ implies a new lower/upper bound for \f$ x_j\f$) */
   SCIP_HASHMAP*         implhash,           /**< hash map from variable to node in implication graph */
   SCIP_Bool**           adjacencymatrix,    /**< adjacencymatrix of conflict graph */
   SCIP_VAR**            totalvars,          /**< problem and SOS1 vars */
   int                   ntotalvars,         /**< number of problem and SOS1 variables*/
   int                   nsos1vars,          /**< number of SOS1 variables */
   int*                  nchgbds,            /**< pointer to store number of changed bounds */
   SCIP_Bool*            implupdate,         /**< pointer to store whether the implication graph has been updated in this function call */
   SCIP_Bool*            cutoff              /**< pointer to store if current nodes LP is infeasible  */
   )
{
   SCIP_CONSHDLR* conshdlrlinear;
   SCIP_CONS** linearconss;
   int nlinearconss;

   SCIP_Bool* implnodes = NULL;     /* implnodes[i] = TRUE if the SOS1 variable corresponding to node i in the implication graph is implied to be nonzero */
   SCIP_Bool* coveredvars = NULL;   /* coveredvars[i] = TRUE if variable with index i is covered by the clique cover */
   int* varindincons = NULL;        /* varindincons[i] = position of SOS1 index i in linear constraint (-1 if x_i is not involved in linear constraint) */

   SCIP_VAR** trafolinvars = NULL;  /* variables of transformed linear constraints without (multi)aggregated variables */
   int ntrafolinvars = 0;
   SCIP_Real* trafolinvals = NULL;
   SCIP_Real* trafoubs = NULL;
   SCIP_Real* trafolbs = NULL;
   SCIP_Real traforhs;
   SCIP_Real trafolhs;

   SCIP_VAR** sos1linvars = NULL;  /* variables that are not contained in linear constraint, but are in conflict with a variable from the linear constraint */
   int nsos1linvars;
   int c;

   assert( scip != NULL );
   assert( conflictgraph != NULL );
   assert( adjacencymatrix != NULL );
   assert( nchgbds != NULL );
   assert( cutoff != NULL );

   *cutoff = FALSE;
   *implupdate = FALSE;

   /* get constraint handler data of linear constraints */
   conshdlrlinear = SCIPfindConshdlr(scip, "linear");
   if ( conshdlrlinear == NULL )
      return SCIP_OKAY;

   /* get linear constraints and number of linear constraints */
   nlinearconss = SCIPconshdlrGetNConss(conshdlrlinear);
   linearconss = SCIPconshdlrGetConss(conshdlrlinear);

   /* allocate buffer arrays */
   SCIP_CALL( SCIPallocBufferArray(scip, &sos1linvars, nsos1vars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &implnodes, nsos1vars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &varindincons, nsos1vars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &coveredvars, ntotalvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &trafoubs, ntotalvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &trafolbs, ntotalvars) );

   /* for every linear constraint and every SOS1 variable */
   for (c = 0; c < nlinearconss + nsos1vars && ! (*cutoff); ++c)
   {
      SCIP_DIGRAPH* conflictgraphlin;
      int** cliquecovers = NULL;           /* clique covers of indices of variables in linear constraint */
      int* cliquecoversizes = NULL;        /* size of each cover */
      int ncliquecovers;
      SCIP_Real* cliquecovervals = NULL;
      int* varincover = NULL;              /* varincover[i] = cover of SOS1 index i */

      int v;
      int i;
      int j;

      /* get transformed linear constraints (without aggregated variables) */
      if ( c < nlinearconss )
      {
         SCIP_VAR** origlinvars;
         int noriglinvars;
         SCIP_Real* origlinvals;
         SCIP_Real origrhs;
         SCIP_Real origlhs;
         SCIP_Real constant;
         int requiredsize;

         /* get data of linear constraint */
         noriglinvars = SCIPgetNVarsLinear(scip, linearconss[c]);
         origlinvars = SCIPgetVarsLinear(scip, linearconss[c]);
         origlinvals = SCIPgetValsLinear(scip, linearconss[c]);
         origrhs = SCIPgetRhsLinear(scip, linearconss[c]);
         origlhs = SCIPgetLhsLinear(scip, linearconss[c]);

         /* allocate buffer arrays */
         SCIP_CALL( SCIPallocBufferArray(scip, &trafolinvars, noriglinvars) );
         SCIP_CALL( SCIPallocBufferArray(scip, &trafolinvals, noriglinvars) );

         /* copy variables and coefficients of linear constraint */
         for (v = 0; v < noriglinvars; ++v)
         {
            trafolinvars[v] = origlinvars[v];
            trafolinvals[v] = origlinvals[v];
         }
         ntrafolinvars = noriglinvars;

         /* transform linear constraint */
         constant = 0.0;
         SCIP_CALL( SCIPgetProbvarLinearSum(scip, trafolinvars, trafolinvals, &ntrafolinvars, noriglinvars, &constant, &requiredsize, TRUE) );
         if( requiredsize > ntrafolinvars )
         {
            SCIP_CALL( SCIPreallocBufferArray(scip, &trafolinvars, requiredsize) );
            SCIP_CALL( SCIPreallocBufferArray(scip, &trafolinvals, requiredsize) );

            SCIP_CALL( SCIPgetProbvarLinearSum(scip, trafolinvars, trafolinvals, &ntrafolinvars, requiredsize, &constant, &requiredsize, TRUE) );
            assert( requiredsize <= ntrafolinvars );
         }
         trafolhs = origlhs - constant;
         traforhs = origrhs - constant;
      }
      else
      {
         SCIP_VAR* var;

         var = nodeGetVarSOS1(conflictgraph, c-nlinearconss);

         if ( SCIPvarGetStatus(var) == SCIP_VARSTATUS_AGGREGATED )
         {
            SCIP_Real constant;

            SCIP_CALL( SCIPallocBufferArray(scip, &trafolinvars, 2) );
            SCIP_CALL( SCIPallocBufferArray(scip, &trafolinvals, 2) );

            constant = SCIPvarGetAggrConstant(var);
            trafolinvars[0] = SCIPvarGetAggrVar(var);
            trafolinvals[0] = SCIPvarGetAggrScalar(var);
            trafolinvars[1] = var;
            trafolinvals[1] = -1.0;
            trafolhs = -constant;
            traforhs = -constant;
            ntrafolinvars = 2;
         }
         else if ( SCIPvarGetStatus(var) == SCIP_VARSTATUS_MULTAGGR )
         {
            SCIP_Real* scalars;
            SCIP_VAR** agrvars;
            SCIP_Real constant;
            int nagrvars;

            nagrvars = SCIPvarGetMultaggrNVars(var);

            SCIP_CALL( SCIPallocBufferArray(scip, &trafolinvars, nagrvars+1) );
            SCIP_CALL( SCIPallocBufferArray(scip, &trafolinvals, nagrvars+1) );

            agrvars = SCIPvarGetMultaggrVars(var);
            scalars = SCIPvarGetMultaggrScalars(var);
            constant = SCIPvarGetMultaggrConstant(var);

            for (v = 0; v < nagrvars; ++v)
            {
               trafolinvars[v] = agrvars[v];
               trafolinvals[v] = scalars[v];
            }
            trafolinvars[nagrvars] = var;
            trafolinvals[nagrvars] = -1.0;
            trafolhs = -constant;
            traforhs = -constant;
            ntrafolinvars = nagrvars + 1;
         }
         else if ( SCIPvarGetStatus(var) == SCIP_VARSTATUS_NEGATED )
         {
            SCIP_VAR* negvar;
            SCIP_Real negcons;

            /* get negation variable and negation offset */
            negvar = SCIPvarGetNegationVar(var);
            negcons = SCIPvarGetNegationConstant(var);

            SCIP_CALL( SCIPallocBufferArray(scip, &trafolinvars, 2) );
            SCIP_CALL( SCIPallocBufferArray(scip, &trafolinvals, 2) );

            trafolinvars[0] = negvar;
            trafolinvars[1] = var;
            trafolinvals[0] = 1.0;
            trafolinvals[1] = 1.0;
            trafolhs = -negcons;
            traforhs = -negcons;
         }
         else
            continue;
      }

      if ( ntrafolinvars == 0 )
      {
         SCIPfreeBufferArray(scip, &trafolinvars);
         SCIPfreeBufferArray(scip, &trafolinvals);
         continue;
      }

      /* compute lower and upper bounds of each term a_i * x_i of transformed constraint */
      for (v = 0; v < ntrafolinvars; ++v)
      {
         SCIP_Real lb = SCIPvarGetLbLocal(trafolinvars[v]);
         SCIP_Real ub = SCIPvarGetUbLocal(trafolinvars[v]);

         if ( trafolinvals[v] < 0.0 )
         {
            SCIP_Real tmp;
            tmp = lb;
            lb = ub;
            ub = tmp;
         }

         if ( SCIPisInfinity(scip, REALABS(lb)) )
            trafolbs[v] = -SCIPinfinity(scip);
         else
            trafolbs[v] = lb * trafolinvals[v];

         if ( SCIPisInfinity(scip, REALABS(ub)) )
            trafoubs[v] = SCIPinfinity(scip);
         else
            trafoubs[v] = ub * trafolinvals[v];
      }

      /* initialization: mark all the SOS1 variables as 'not a member of the linear constraint' */
      for (v = 0; v < nsos1vars; ++v)
         varindincons[v] = -1;

      /* save position of SOS1 variables in linear constraint */
      for (v = 0; v < ntrafolinvars; ++v)
      {
         int node;

         node = varGetNodeSOS1(conshdlr, trafolinvars[v]);

         if ( node >= 0 )
            varindincons[node] = v;
      }

      /* create conflict graph of linear constraint */
      SCIP_CALL( SCIPdigraphCreate(&conflictgraphlin, ntrafolinvars) );
      SCIP_CALL( SCIPgenConflictgraphLinearCons(conshdlr, conflictgraphlin, conflictgraph, trafolinvars, ntrafolinvars, varindincons) );

      /* mark all the variables as 'not covered by some clique cover' */
      for (i = 0; i < ntrafolinvars; ++i)
         coveredvars[i] = FALSE;

      /* allocate buffer array */
      SCIP_CALL( SCIPallocBufferArray(scip, &cliquecovervals, ntrafolinvars) );
      SCIP_CALL( SCIPallocBufferArray(scip, &cliquecoversizes, ntrafolinvars) );
      SCIP_CALL( SCIPallocBufferArray(scip, &cliquecovers, ntrafolinvars) );

      /* compute distinct cliques that cover all the variables of the linear constraint */
      ncliquecovers = 0;
      for (v = 0; v < ntrafolinvars; ++v)
      {
         /* if variable is not already covered by an already known clique cover */
         if ( ! coveredvars[v] )
         {
            SCIP_CALL( SCIPallocBufferArray(scip, &(cliquecovers[ncliquecovers]), ntrafolinvars) );
            SCIP_CALL( SCIPcomputeVarsCoverSOS1(scip, conflictgraph, conflictgraphlin, trafolinvars, coveredvars, cliquecovers[ncliquecovers], &(cliquecoversizes[ncliquecovers]), v, FALSE) );
            ++ncliquecovers;
         }
      }

      /* free conflictgraph */
      SCIPdigraphFree(&conflictgraphlin);

      /* compute variables that are not contained in transformed linear constraint, but are in conflict with a variable from the transformed linear constraint */
      nsos1linvars = 0;
      for (v = 0; v < ntrafolinvars; ++v)
      {
         int nodev;

         nodev = varGetNodeSOS1(conshdlr, trafolinvars[v]);

         /* if variable is an SOS1 variable */
         if ( nodev >= 0 )
         {
            int succnode;
            int nsucc;
            int* succ;
            int s;

            succ = SCIPdigraphGetSuccessors(conflictgraph, nodev);
            nsucc = SCIPdigraphGetNSuccessors(conflictgraph, nodev);

            for (s = 0; s < nsucc; ++s)
            {
               succnode = succ[s];

               /* if variable is not a member of linear constraint and not already listed in the array sos1linvars */
               if ( varindincons[succnode] == -1 )
               {
                  sos1linvars[nsos1linvars] = nodeGetVarSOS1(conflictgraph, succnode);
                  varindincons[succnode] = -2; /* mark variable as listed in array sos1linvars */
                  ++nsos1linvars;
               }
            }
         }
      }


      /* try to tighten lower bounds */

      /* sort each cliquecover array in ascending order of the lower bounds of a_i * x_i; fill vector varincover */
      SCIP_CALL( SCIPallocBufferArray(scip, &varincover, ntrafolinvars) );
      for (i = 0; i < ncliquecovers; ++i)
      {
         for (j = 0; j < cliquecoversizes[i]; ++j)
         {
            int ind = cliquecovers[i][j];

            varincover[ind] = i;
            cliquecovervals[j] = trafoubs[ind];
         }
         SCIPsortDownRealInt(cliquecovervals, cliquecovers[i], cliquecoversizes[i]);
      }

      /* for every variable in transformed constraint: try lower bound tightening */
      for (v = 0; v < ntrafolinvars + nsos1linvars; ++v)
      {
         SCIP_Real newboundnonzero; /* new bound of a_v*x_v if we assume that x_v != 0 */
         SCIP_Real newboundnores;   /* new bound of a_v*x_v if we assume that x_v = 0 is possible */
         SCIP_Real newbound;        /* resulting new bound of x_v */
         SCIP_VAR* var;
         SCIP_Real trafoubv;
         SCIP_Real linval;
         SCIP_Real ub;
         SCIP_Real lb;
         SCIP_Bool tightened;
         SCIP_Bool infeasible;
         SCIP_Bool inftynores = FALSE;
         SCIP_Bool update;
         int ninftynonzero = 0;
         int nodev;
         int w;

         if ( v < ntrafolinvars )
         {
            var = trafolinvars[v];
            trafoubv = trafoubs[v];
         }
         else
         {
            assert( v-ntrafolinvars >= 0 );
            var = sos1linvars[v-ntrafolinvars];/*lint !e679*/
            trafoubv = 0.0;
         }

         ub = SCIPvarGetUbLocal(var);
         lb = SCIPvarGetLbLocal(var);

         if ( SCIPisInfinity(scip, -trafolhs) || SCIPisZero(scip, ub - lb) )
            continue;

         newboundnonzero = trafolhs;
         newboundnores = trafolhs;
         nodev = varGetNodeSOS1(conshdlr, var); /* possibly -1 if var is not involved in an SOS1 constraint */
         assert( nodev < nsos1vars );

         /* determine incidence vector of implication variables */
         for (w = 0; w < nsos1vars; ++w)
            implnodes[w] = FALSE;
         SCIP_CALL( SCIPgetSOS1Implications(scip, conshdlr, totalvars, implgraph, implhash, implnodes, (int) (size_t) SCIPhashmapGetImage(implhash, var)) );

         /* compute new bound */
         for (i = 0; i < ncliquecovers; ++i)
         {
            int indcliq;
            int nodecliq;

            assert( cliquecoversizes[i] > 0 );

            indcliq = cliquecovers[i][0];
            assert( 0 <= indcliq && indcliq < ntrafolinvars );

            /* determine maximum without index v (note that trafoub is sorted non-increasingly) */
            if ( v != indcliq )
            {
               if ( SCIPisInfinity(scip, trafoubs[indcliq]) )
                  inftynores = TRUE;
               else
                  newboundnores -= trafoubs[indcliq];
            }
            else if ( cliquecoversizes[i] > 1 )
            {
               assert( 0 <= cliquecovers[i][1] && cliquecovers[i][1] < ntrafolinvars );
               if ( SCIPisInfinity(scip, trafoubs[cliquecovers[i][1]]) )
                  inftynores = TRUE;
               else
                  newboundnores -= trafoubs[cliquecovers[i][1]];/*lint --e{679}*/
            }

            for (j = 0; j < cliquecoversizes[i]; ++j)
            {
               indcliq = cliquecovers[i][j];
               assert( 0 <= indcliq && indcliq < ntrafolinvars );

               nodecliq = varGetNodeSOS1(conshdlr, trafolinvars[indcliq]); /* possibly -1 if variable is not involved in an SOS1 constraint */
               assert( nodecliq < nsos1vars );

               if ( v != indcliq )
               {
                  /* if nodev or nodecliq are not a member of an SOS1 constraint or the variable corresponding to nodecliq is not implied to be zero if x_v != 0  */
                  if ( nodev < 0 || nodecliq < 0 || (! isConnectedSOS1(adjacencymatrix, NULL, nodev, nodecliq) && ! isImpliedZero(conflictgraph, implnodes, nodecliq) ) )
                  {
                     if ( SCIPisInfinity(scip, trafoubs[indcliq]) )
                        ++ninftynonzero;
                     else
                        newboundnonzero -= trafoubs[indcliq];
                     break;
                  }
               }
            }
         }
         assert( ninftynonzero == 0 || inftynores );

         /* if computed upper bound is not infinity and variable is contained in linear constraint */
         if ( ninftynonzero == 0 && v < ntrafolinvars )
         {
            linval = trafolinvals[v];
            /* compute new bound */
            assert( ! SCIPisFeasZero(scip, linval) );
            if ( SCIPisFeasPositive(scip, newboundnores) && ! inftynores )
               newbound = newboundnonzero;
            else
               newbound = MIN(0, newboundnonzero);
            newbound /= linval;

            /* check if new bound is tighter than the old one or problem is infeasible */
            if ( SCIPisFeasPositive(scip, linval) && SCIPisFeasLT(scip, lb, newbound) )
            {
               if ( SCIPisFeasLT(scip, ub, newbound) )
               {
                  *cutoff = TRUE;
                  break;
               }

               if ( SCIPvarIsIntegral(var) )
                  newbound = SCIPceil(scip, newbound);

               SCIP_CALL( SCIPtightenVarLb(scip, var, newbound, FALSE, &infeasible, &tightened) );
               assert( ! infeasible );

               if ( tightened )
               {
                  SCIPdebugMessage("changed lower bound of variable %s from %f to %f \n", SCIPvarGetName(var), lb, newbound);
                  ++(*nchgbds);
               }
            }
            else if ( SCIPisFeasNegative(scip, linval) && SCIPisFeasGT(scip, ub, newbound) )
            {
               /* if assumption a_i * x_i != 0 was not correct */
               if ( SCIPisFeasGT(scip, SCIPvarGetLbLocal(var), newbound) )
               {
                  *cutoff = TRUE;
                  break;
               }

               if ( SCIPvarIsIntegral(var) )
                  newbound = SCIPfloor(scip, newbound);

               SCIP_CALL( SCIPtightenVarUb(scip, var, newbound, FALSE, &infeasible, &tightened) );
               assert( ! infeasible );

               if ( tightened )
               {
                  SCIPdebugMessage("changed upper bound of variable %s from %f to %f \n", SCIPvarGetName(var), ub, newbound);
                  ++(*nchgbds);
               }
            }
         }

         /* update implication graph if possible */
         SCIP_CALL( updateImplicationGraphSOS1(scip, conshdlr, conflictgraph, adjacencymatrix, implgraph, implhash, implnodes, totalvars, cliquecovers, cliquecoversizes, varincover,
               trafolinvars, trafolinvals, ntrafolinvars, trafoubs, var, trafoubv, newboundnonzero, ninftynonzero, TRUE, nchgbds, &update) );
         if ( update )
            *implupdate = TRUE;
      }

      if ( *cutoff == TRUE )
      {
         /* free memory */
         SCIPfreeBufferArrayNull(scip, &varincover);
         for (j = ncliquecovers-1; j >= 0; --j)
            SCIPfreeBufferArrayNull(scip, &cliquecovers[j]);
         SCIPfreeBufferArrayNull(scip, &cliquecovers);
         SCIPfreeBufferArrayNull(scip, &cliquecoversizes);
         SCIPfreeBufferArrayNull(scip, &cliquecovervals);
         SCIPfreeBufferArrayNull(scip, &trafolinvals);
         SCIPfreeBufferArrayNull(scip, &trafolinvars);
         break;
      }


      /* try to tighten upper bounds */

      /* sort each cliquecover array in ascending order of the lower bounds of a_i * x_i; fill vector varincover */
      for (i = 0; i < ncliquecovers; ++i)
      {
         for (j = 0; j < cliquecoversizes[i]; ++j)
         {
            int ind = cliquecovers[i][j];

            varincover[ind] = i;
            cliquecovervals[j] = trafolbs[ind];
         }
         SCIPsortRealInt(cliquecovervals, cliquecovers[i], cliquecoversizes[i]);
      }

      /* for every variable that is in transformed constraint or every variable that is in conflict with some variable from trans. cons.:
         try upper bound tightening */
      for (v = 0; v < ntrafolinvars + nsos1linvars; ++v)
      {
         SCIP_Real newboundnonzero; /* new bound of a_v*x_v if we assume that x_v != 0 */
         SCIP_Real newboundnores;   /* new bound of a_v*x_v if there are no restrictions */
         SCIP_Real newbound;        /* resulting new bound of x_v */
         SCIP_VAR* var;
         SCIP_Real linval;
         SCIP_Real trafolbv;
         SCIP_Real lb;
         SCIP_Real ub;
         SCIP_Bool tightened;
         SCIP_Bool infeasible;
         SCIP_Bool inftynores = FALSE;
         SCIP_Bool update;
         int ninftynonzero = 0;
         int nodev;
         int w;

         if ( v < ntrafolinvars )
         {
            var = trafolinvars[v];
            trafolbv = trafolbs[v];
         }
         else
         {
            assert( v-ntrafolinvars >= 0 );
            var = sos1linvars[v-ntrafolinvars];/*lint !e679*/
            trafolbv = 0.0; /* since variable is not a member of linear constraint */
         }
         lb = SCIPvarGetLbLocal(var);
         ub = SCIPvarGetUbLocal(var);
         if ( SCIPisInfinity(scip, traforhs) || SCIPisEQ(scip, lb, ub) )
            continue;

         newboundnonzero = traforhs;
         newboundnores = traforhs;
         nodev = varGetNodeSOS1(conshdlr, var); /* possibly -1 if var is not involved in an SOS1 constraint */
         assert( nodev < nsos1vars );

         /* determine incidence vector of implication variables (i.e., which SOS1 variables are nonzero if x_v is nonzero) */
         for (w = 0; w < nsos1vars; ++w)
            implnodes[w] = FALSE;
         SCIP_CALL( SCIPgetSOS1Implications(scip, conshdlr, totalvars, implgraph, implhash, implnodes, (int) (size_t) SCIPhashmapGetImage(implhash, var)) );

         /* compute new bound */
         for (i = 0; i < ncliquecovers; ++i)
         {
            int indcliq;
            int nodecliq;

            assert( cliquecoversizes[i] > 0 );

            indcliq = cliquecovers[i][0];
            assert( 0 <= indcliq && indcliq < ntrafolinvars );

            /* determine maximum without index v (note that trafolbs is sorted non-increasingly) */
            if ( v != indcliq )
            {
               /* if bound would be infinity */
               if ( SCIPisInfinity(scip, -trafolbs[indcliq]) )
                  inftynores = TRUE;
               else
                  newboundnores -= trafolbs[indcliq];
            }
            else if ( cliquecoversizes[i] > 1 )
            {
               assert( 0 <= cliquecovers[i][1] && cliquecovers[i][1] < ntrafolinvars );
               if ( SCIPisInfinity(scip, -trafolbs[cliquecovers[i][1]]) )
                  inftynores = TRUE;
               else
                  newboundnores -= trafolbs[cliquecovers[i][1]]; /*lint --e{679}*/
            }

            /* determine maximum without index v and if x_v is nonzero (note that trafolbs is sorted non-increasingly) */
            for (j = 0; j < cliquecoversizes[i]; ++j)
            {
               indcliq = cliquecovers[i][j];
               assert( 0 <= indcliq && indcliq < ntrafolinvars );

               nodecliq = varGetNodeSOS1(conshdlr, trafolinvars[indcliq]); /* possibly -1 if variable is not involved in an SOS1 constraint */
               assert( nodecliq < nsos1vars );

               if ( v != indcliq )
               {
                  /* if nodev or nodecliq are not a member of an SOS1 constraint or the variable corresponding to nodecliq is not implied to be zero if x_v != 0  */
                  if ( nodev < 0 || nodecliq < 0 || (! isConnectedSOS1(adjacencymatrix, NULL, nodev, nodecliq) && ! isImpliedZero(conflictgraph, implnodes, nodecliq) ) )
                  {
                     /* if bound would be infinity */
                     if ( SCIPisInfinity(scip, -trafolbs[indcliq]) )
                        ++ninftynonzero;
                     else
                        newboundnonzero -= trafolbs[indcliq];
                     break;
                  }
               }
            }
         }
         assert( ninftynonzero == 0 || inftynores );


         /* if computed upper bound is not infinity and variable is contained in linear constraint */
         if ( ninftynonzero == 0 && v < ntrafolinvars )
         {
            linval = trafolinvals[v];

            /* compute new bound */
            assert( ! SCIPisFeasZero(scip, linval) );
            if ( SCIPisFeasNegative(scip, newboundnores) && ! inftynores )
               newbound = newboundnonzero;
            else
               newbound = MAX(0, newboundnonzero);
            newbound /= linval;

            /* check if new bound is tighter than the old one or problem is infeasible */
            if ( SCIPisFeasPositive(scip, linval) && SCIPisFeasGT(scip, ub, newbound) )
            {
               /* if new upper bound is smaller than the lower bound, we are infeasible */
               if ( SCIPisFeasGT(scip, lb, newbound) )
               {
                  *cutoff = TRUE;
                  break;
               }

               if ( SCIPvarIsIntegral(var) )
                  newbound = SCIPfloor(scip, newbound);

               SCIP_CALL( SCIPtightenVarUb(scip, var, newbound, FALSE, &infeasible, &tightened) );
               assert( ! infeasible );

               if ( tightened )
               {
                  SCIPdebugMessage("changed upper bound of variable %s from %f to %f \n", SCIPvarGetName(var), ub, newbound);
                  ++(*nchgbds);
               }
            }
            else if ( SCIPisFeasNegative(scip, linval) && SCIPisFeasLT(scip, lb, newbound) )
            {
               /* if assumption a_i * x_i != 0 was not correct */
               if ( SCIPisFeasLT(scip, ub, newbound) )
               {
                  *cutoff = TRUE;
                  break;
               }

               if ( SCIPvarIsIntegral(var) )
                  newbound = SCIPceil(scip, newbound);

               SCIP_CALL( SCIPtightenVarLb(scip, var, newbound, FALSE, &infeasible, &tightened) );
               assert( ! infeasible );

               if ( tightened )
               {
                  SCIPdebugMessage("changed lower bound of variable %s from %f to %f \n", SCIPvarGetName(var), lb, newbound);
                  ++(*nchgbds);
               }
            }
         }

         /* update implication graph if possible */
         SCIP_CALL( updateImplicationGraphSOS1(scip, conshdlr, conflictgraph, adjacencymatrix, implgraph, implhash, implnodes, totalvars, cliquecovers, cliquecoversizes, varincover,
               trafolinvars, trafolinvals, ntrafolinvars, trafolbs, var, trafolbv, newboundnonzero, ninftynonzero, FALSE, nchgbds, &update) );
         if ( update )
            *implupdate = TRUE;
      }

      /* free memory */
      SCIPfreeBufferArrayNull(scip, &varincover);
      for (j = ncliquecovers-1; j >= 0; --j)
         SCIPfreeBufferArrayNull(scip, &cliquecovers[j]);
      SCIPfreeBufferArrayNull(scip, &cliquecovers);
      SCIPfreeBufferArrayNull(scip, &cliquecoversizes);
      SCIPfreeBufferArrayNull(scip, &cliquecovervals);
      SCIPfreeBufferArrayNull(scip, &trafolinvals);
      SCIPfreeBufferArrayNull(scip, &trafolinvars);

      if ( *cutoff == TRUE )
         break;
   } /* end for every linear constraint */

   /* free buffer arrays */
   SCIPfreeBufferArrayNull(scip, &sos1linvars);
   SCIPfreeBufferArrayNull(scip, &trafolbs);
   SCIPfreeBufferArrayNull(scip, &trafoubs);
   SCIPfreeBufferArrayNull(scip, &coveredvars);
   SCIPfreeBufferArrayNull(scip, &varindincons);
   SCIPfreeBufferArrayNull(scip, &implnodes);

   return SCIP_OKAY;
}


/* ----------------------------- propagation -------------------------------------*/

/** propagate variables */
static
SCIP_RETCODE propSOS1(
   SCIP*                 scip,               /**< SCIP pointer */
   SCIP_CONS*            cons,               /**< constraint */
   SCIP_CONSDATA*        consdata,           /**< constraint data */
   SCIP_Bool*            cutoff,             /**< whether a cutoff happened */
   int*                  ngen                /**< pointer to incremental counter for domain changes */
   )
{
   assert( scip != NULL );
   assert( cons != NULL );
   assert( consdata != NULL );
   assert( cutoff != NULL );
   assert( ngen != NULL );

   *cutoff = FALSE;

   /* if more than one variable is fixed to be nonzero */
   if ( consdata->nfixednonzeros > 1 )
   {
      SCIPdebugMessage("the node is infeasible, more than 1 variable is fixed to be nonzero.\n");
      SCIP_CALL( SCIPresetConsAge(scip, cons) );
      *cutoff = TRUE;
      return SCIP_OKAY;
   }

   /* if exactly one variable is fixed to be nonzero */
   if ( consdata->nfixednonzeros == 1 )
   {
      SCIP_VAR** vars;
      SCIP_Bool infeasible;
      SCIP_Bool tightened;
      SCIP_Bool success;
      SCIP_Bool allVarFixed;
      int firstFixedNonzero;
      int ngenold;
      int nvars;
      int j;

      firstFixedNonzero = -1;
      nvars = consdata->nvars;
      vars = consdata->vars;
      assert( vars != NULL );
      ngenold = *ngen;

      /* search nonzero variable - is needed for propinfo */
      for (j = 0; j < nvars; ++j)
      {
         if ( SCIPisFeasPositive(scip, SCIPvarGetLbLocal(vars[j])) || SCIPisFeasNegative(scip, SCIPvarGetUbLocal(vars[j])) )
         {
            firstFixedNonzero = j;
            break;
         }
      }
      assert( firstFixedNonzero >= 0 );

      SCIPdebugMessage("variable <%s> is fixed nonzero, fixing other variables to 0.\n", SCIPvarGetName(vars[firstFixedNonzero]));

      /* fix variables before firstFixedNonzero to 0 */
      allVarFixed = TRUE;
      for (j = 0; j < firstFixedNonzero; ++j)
      {
         /* fix variable */
         SCIP_CALL( inferVariableZero(scip, vars[j], cons, firstFixedNonzero, &infeasible, &tightened, &success) );
         assert( ! infeasible );
         allVarFixed = allVarFixed && success;
         if ( tightened )
            ++(*ngen);
      }

      /* fix variables after firstFixedNonzero to 0 */
      for (j = firstFixedNonzero+1; j < nvars; ++j)
      {
         /* fix variable */
         SCIP_CALL( inferVariableZero(scip, vars[j], cons, firstFixedNonzero, &infeasible, &tightened, &success) );
         assert( ! infeasible ); /* there should be no variables after firstFixedNonzero that are fixed to be nonzero */
         allVarFixed = allVarFixed && success;
         if ( tightened )
            ++(*ngen);
      }

      /* reset constraint age counter */
      if ( *ngen > ngenold )
      {
         SCIP_CALL( SCIPresetConsAge(scip, cons) );
      }

      /* delete constraint locally */
      if ( allVarFixed )
      {
         assert( !SCIPconsIsModifiable(cons) );
         SCIP_CALL( SCIPdelConsLocal(scip, cons) );
      }
   }

   return SCIP_OKAY;
}

/* ----------------------------- branching -------------------------------------*/

/** enforcement method
 *
 *  We check whether the current solution is feasible, i.e., contains at most one nonzero
 *  variable. If not, we branch along the lines indicated by Beale and Tomlin:
 *
 *  We first compute \f$W = \sum_{j=1}^n |x_i|\f$ and \f$w = \sum_{j=1}^n j\, |x_i|\f$. Then we
 *  search for the index \f$k\f$ that satisfies
 *  \f[
 *        k \leq \frac{w}{W} < k+1.
 *  \f]
 *  The branches are then
 *  \f[
 *        x_1 = 0, \ldots, x_k = 0 \qquad \mbox{and}\qquad x_{k+1} = 0, \ldots, x_n = 0.
 *  \f]
 *
 *  If the constraint contains two variables, the branching of course simplifies.
 *
 *  Depending on the parameters (@c branchnonzeros, @c branchweight) there are three ways to choose
 *  the branching constraint.
 *
 *  <TABLE>
 *  <TR><TD>@c branchnonzeros</TD><TD>@c branchweight</TD><TD>constraint chosen</TD></TR>
 *  <TR><TD>@c true          </TD><TD> ?             </TD><TD>most number of nonzeros</TD></TR>
 *  <TR><TD>@c false         </TD><TD> @c true       </TD><TD>maximal weight corresponding to nonzero variable</TD></TR>
 *  <TR><TD>@c false         </TD><TD> @c true       </TD><TD>largest sum of variable values</TD></TR>
 *  </TABLE>
 *
 *  @c branchnonzeros = @c false, @c branchweight = @c true allows the user to specify an order for
 *  the branching importance of the constraints (setting the weights accordingly).
 *
 *  Constraint branching can also be turned off using parameter @c branchsos.
 */
static
SCIP_RETCODE enforceSOS1(
   SCIP*                 scip,               /**< SCIP pointer */
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler */
   int                   nconss,             /**< number of constraints */
   SCIP_CONS**           conss,              /**< indicator constraints */
   SCIP_RESULT*          result              /**< result */
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONSDATA* consdata;
   SCIP_NODE* node1;
   SCIP_NODE* node2;
   SCIP_CONS* branchCons;
   SCIP_Real maxWeight;
   SCIP_VAR** vars;
   int nvars;
   int c;

   assert( scip != NULL );
   assert( conshdlr != NULL );
   assert( conss != NULL );
   assert( result != NULL );

   maxWeight = -SCIP_REAL_MAX;
   branchCons = NULL;

   SCIPdebugMessage("Enforcing SOS1 constraints <%s>.\n", SCIPconshdlrGetName(conshdlr) );
   *result = SCIP_FEASIBLE;

   /* get constraint handler data */
   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert( conshdlrdata != NULL );

   /* check each constraint */
   for (c = 0; c < nconss; ++c)
   {
      SCIP_CONS* cons;
      SCIP_Bool cutoff;
      SCIP_Real weight;
      int ngen;
      int cnt;
      int j;

      cons = conss[c];
      assert( cons != NULL );
      consdata = SCIPconsGetData(cons);
      assert( consdata != NULL );

      ngen = 0;
      cnt = 0;
      nvars = consdata->nvars;
      vars = consdata->vars;

      /* do nothing if there are not enough variables - this is usually eliminated by preprocessing */
      if ( nvars < 2 )
         continue;

      /* first perform propagation (it might happen that standard propagation is turned off) */
      SCIP_CALL( propSOS1(scip, cons, consdata, &cutoff, &ngen) );
      SCIPdebugMessage("propagating <%s> in enforcing (cutoff: %u, domain reductions: %d).\n", SCIPconsGetName(cons), cutoff, ngen);
      if ( cutoff )
      {
         *result = SCIP_CUTOFF;
         return SCIP_OKAY;
      }
      if ( ngen > 0 )
      {
         *result = SCIP_REDUCEDDOM;
         return SCIP_OKAY;
      }
      assert( ngen == 0 );

      /* check constraint */
      weight = 0.0;
      for (j = 0; j < nvars; ++j)
      {
         SCIP_Real val = REALABS(SCIPgetSolVal(scip, NULL, vars[j]));

         if ( ! SCIPisFeasZero(scip, val) )
         {
            if ( conshdlrdata->branchnonzeros )
               weight += 1.0;
            else
            {
               if ( conshdlrdata->branchweight )
               {
                  /* choose maximum nonzero-variable weight */
                  if ( consdata->weights[j] > weight )
                     weight = consdata->weights[j];
               }
               else
                  weight += val;
            }
            ++cnt;
         }
      }
      /* if constraint is violated */
      if ( cnt > 1 && weight > maxWeight )
      {
         maxWeight = weight;
         branchCons = cons;
      }
   }

   /* if all constraints are feasible */
   if ( branchCons == NULL )
   {
      SCIPdebugMessage("All SOS1 constraints are feasible.\n");
      return SCIP_OKAY;
   }

   /* if we should leave branching decision to branching rules */
   if ( ! conshdlrdata->branchsos )
   {
      *result = SCIP_INFEASIBLE;
      return SCIP_OKAY;
   }

   /* otherwise create branches */
   SCIPdebugMessage("Branching on constraint <%s> (weight: %f).\n", SCIPconsGetName(branchCons), maxWeight);
   consdata = SCIPconsGetData(branchCons);
   assert( consdata != NULL );
   nvars = consdata->nvars;
   vars = consdata->vars;

   if ( nvars == 2 )
   {
      SCIP_Bool infeasible;

      /* constraint is infeasible: */
      assert( ! SCIPisFeasZero(scip, SCIPgetSolVal(scip, NULL, vars[0])) && ! SCIPisFeasZero(scip, SCIPgetSolVal(scip, NULL, vars[1])) );

      /* create branches */
      SCIPdebugMessage("Creating two branches.\n");

      SCIP_CALL( SCIPcreateChild(scip, &node1, SCIPcalcNodeselPriority(scip, vars[0], SCIP_BRANCHDIR_DOWNWARDS, 0.0), SCIPcalcChildEstimate(scip, vars[0], 0.0) ) );
      SCIP_CALL( fixVariableZeroNode(scip, vars[0], node1, &infeasible) );
      assert( ! infeasible );

      SCIP_CALL( SCIPcreateChild(scip, &node2, SCIPcalcNodeselPriority(scip, vars[1], SCIP_BRANCHDIR_DOWNWARDS, 0.0), SCIPcalcChildEstimate(scip, vars[1], 0.0) ) );
      SCIP_CALL( fixVariableZeroNode(scip, vars[1], node2, &infeasible) );
      assert( ! infeasible );
   }
   else
   {
      SCIP_Bool infeasible;
      SCIP_Real weight1;
      SCIP_Real weight2;
      SCIP_Real nodeselest;
      SCIP_Real objest;
      SCIP_Real w;
      int j;
      int ind;
      int cnt;

      cnt = 0;

      weight1 = 0.0;
      weight2 = 0.0;

      /* compute weight */
      for (j = 0; j < nvars; ++j)
      {
         SCIP_Real val = REALABS(SCIPgetSolVal(scip, NULL, vars[j]));
         weight1 += val * (SCIP_Real) j;
         weight2 += val;

         if ( ! SCIPisFeasZero(scip, val) )
            ++cnt;
      }

      assert( cnt >= 2 );
      assert( !SCIPisFeasZero(scip, weight2) );
      w = weight1/weight2;  /*lint !e795*/

      ind = (int) SCIPfloor(scip, w);
      assert( 0 <= ind && ind < nvars-1 );

      /* branch on variable ind: either all variables up to ind or all variables after ind are zero */
      SCIPdebugMessage("Branching on variable <%s>.\n", SCIPvarGetName(vars[ind]));

      /* calculate node selection and objective estimate for node 1 */
      nodeselest = 0.0;
      objest = 0.0;
      for (j = 0; j <= ind; ++j)
      {
         nodeselest += SCIPcalcNodeselPriority(scip, vars[j], SCIP_BRANCHDIR_DOWNWARDS, 0.0);
         objest += SCIPcalcChildEstimate(scip, vars[j], 0.0);
      }
      /* take the average of the individual estimates */
      objest = objest/((SCIP_Real) ind + 1.0);

      /* create node 1 */
      SCIP_CALL( SCIPcreateChild(scip, &node1, nodeselest, objest) );
      for (j = 0; j <= ind; ++j)
      {
         SCIP_CALL( fixVariableZeroNode(scip, vars[j], node1, &infeasible) );
         assert( ! infeasible );
      }

      /* calculate node selection and objective estimate for node 1 */
      nodeselest = 0.0;
      objest = 0.0;
      for (j = ind+1; j < nvars; ++j)
      {
         nodeselest += SCIPcalcNodeselPriority(scip, vars[j], SCIP_BRANCHDIR_DOWNWARDS, 0.0);
         objest += SCIPcalcChildEstimate(scip, vars[j], 0.0);
      }
      /* take the average of the individual estimates */
      objest = objest/((SCIP_Real) (nvars - ind - 1));

      /* create node 2 */
      SCIP_CALL( SCIPcreateChild(scip, &node2, nodeselest, objest) );
      for (j = ind+1; j < nvars; ++j)
      {
         SCIP_CALL( fixVariableZeroNode(scip, vars[j], node2, &infeasible) );
         assert( ! infeasible );
      }
   }
   SCIP_CALL( SCIPresetConsAge(scip, branchCons) );
   *result = SCIP_BRANCHED;

   return SCIP_OKAY;
}


/* ----------------------------- separation ------------------------------------*/

/* initialitze tclique graph and create clique data */
static
SCIP_RETCODE initTCliquegraph(
   SCIP*                 scip,               /**< SCIP pointer */
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler */
   SCIP_CONSHDLRDATA*    conshdlrdata,       /**< constraint handler data */
   SCIP_DIGRAPH*         conflictgraph,      /**< conflict graph */
   int                   nsos1vars,          /**< number of SOS1 variables */
   SCIP_SOL*             sol,                /**< LP solution to be separated (or NULL) */
   SCIP_Real             scaleval            /**< factor for scaling weights */
   )
{
   TCLIQUE_DATA* tcliquedata;
   int j;

   /* try to generate bound cuts */
   if ( ! tcliqueCreate(&conshdlrdata->tcliquegraph) )
      return SCIP_NOMEMORY;

   /* add nodes */
   for (j = 0; j < nsos1vars; ++j)
   {
      if ( ! tcliqueAddNode(conshdlrdata->tcliquegraph, j, 0 ) )
         return SCIP_NOMEMORY;
   }

   /* add edges */
   for (j = 0; j < nsos1vars; ++j)
   {
      int* succ;
      int nsucc;
      int succnode;
      int i;

      nsucc = SCIPdigraphGetNSuccessors(conflictgraph, j);
      succ = SCIPdigraphGetSuccessors(conflictgraph, j);

      for (i = 0; i < nsucc; ++i)
      {
         succnode = succ[i];

         if ( succnode > j && SCIPvarIsActive(nodeGetVarSOS1(conflictgraph, succnode)) )
         {
            if ( ! tcliqueAddEdge(conshdlrdata->tcliquegraph, j, succnode) )
               return SCIP_NOMEMORY;
         }
      }
   }
   if ( ! tcliqueFlush(conshdlrdata->tcliquegraph) )
      return SCIP_NOMEMORY;


   /* allocate clique data */
   SCIP_CALL( SCIPallocMemory(scip, &conshdlrdata->tcliquedata) );
   tcliquedata = conshdlrdata->tcliquedata;

   /* initialize clique data */
   tcliquedata->scip = scip;
   tcliquedata->sol = sol;
   tcliquedata->conshdlr = conshdlr;
   tcliquedata->conflictgraph = conflictgraph;
   tcliquedata->scaleval = scaleval;
   tcliquedata->ncuts = 0;
   tcliquedata->nboundcuts = conshdlrdata->nboundcuts;
   tcliquedata->strthenboundcuts = conshdlrdata->strthenboundcuts;
   tcliquedata->maxboundcuts = conshdlrdata->maxboundcutsroot;

   return SCIP_OKAY;
}


/* update weights of tclique graph */
static
SCIP_RETCODE updateWeightsTCliquegraph(
   SCIP*                 scip,               /**< SCIP pointer */
   SCIP_CONSHDLRDATA*    conshdlrdata,       /**< constraint handler data */
   TCLIQUE_DATA*         tcliquedata,        /**< tclique data */
   SCIP_DIGRAPH*         conflictgraph,      /**< conflict graph */
   SCIP_SOL*             sol,                /**< LP solution to be separated (or NULL) */
   int                   nsos1vars           /**< number of SOS1 variables */
   )
{
   SCIP_Real scaleval;
   int j;

   scaleval = tcliquedata->scaleval;

   for (j = 0; j < nsos1vars; ++j)
   {
      SCIP_Real solval;
      SCIP_Real bound;
      SCIP_VAR* var;

      var = nodeGetVarSOS1(conflictgraph, j);
      solval = SCIPgetSolVal(scip, sol, var);

      if ( SCIPisFeasPositive(scip, solval) )
      {
         if ( conshdlrdata->strthenboundcuts )
            bound = REALABS( SCIPnodeGetSolvalVarboundUbSOS1(scip, conflictgraph, sol, j) );
         else
            bound = REALABS( SCIPvarGetUbLocal(var) );
      }
      else if ( SCIPisFeasNegative(scip, solval) )
      {
         if ( conshdlrdata->strthenboundcuts )
            bound = REALABS( SCIPnodeGetSolvalVarboundLbSOS1(scip, conflictgraph, sol, j) );
         else
            bound = REALABS( SCIPvarGetLbLocal(var) );
      }
      else
         bound = 0.0;

      solval = REALABS( solval );

      if ( ! SCIPisFeasZero(scip, bound) && ! SCIPisInfinity(scip, bound) )
      {
         SCIP_Real nodeweight = REALABS( solval/bound ) * scaleval;
         tcliqueChangeWeight(conshdlrdata->tcliquegraph, j, (int)nodeweight);
      }
      else
      {
         tcliqueChangeWeight(conshdlrdata->tcliquegraph, j, 0);
      }
   }

   return SCIP_OKAY;
}


/* adds bound cut(s) to separation storage */
static
SCIP_RETCODE addBoundCutSepa(
   SCIP*                 scip,               /**< SCIP pointer */
   TCLIQUE_DATA*         tcliquedata,        /**< clique data */
   SCIP_ROW*             rowlb,              /**< row for lower bounds (or NULL) */
   SCIP_ROW*             rowub,              /**< row for upper bounds (or NULL) */
   SCIP_Bool*            success             /**< pointer to store if bound cut was added */
   )
{
   assert( scip != NULL );
   assert( tcliquedata != NULL );
   assert( success != NULL);

   *success = FALSE;

   /* add cut for lower bounds */
   if ( rowlb != NULL )
   {
      if ( ! SCIProwIsInLP(rowlb) && SCIPisCutEfficacious(scip, NULL, rowlb) )
      {
         SCIP_Bool infeasible;

         SCIP_CALL( SCIPaddCut(scip, NULL, rowlb, FALSE, &infeasible) );
         assert( ! infeasible );
         SCIPdebug( SCIP_CALL( SCIPprintRow(scip, rowlb, NULL) ) );
         ++tcliquedata->nboundcuts;
         ++tcliquedata->ncuts;
         *success = TRUE;
      }
      SCIP_CALL( SCIPreleaseRow(scip, &rowlb) );
   }

   /* add cut for upper bounds */
   if ( rowub != NULL )
   {
      if ( ! SCIProwIsInLP(rowub) && SCIPisCutEfficacious(scip, NULL, rowub) )
      {
         SCIP_Bool infeasible;

         SCIP_CALL( SCIPaddCut(scip, NULL, rowub, FALSE, &infeasible) );
         assert( ! infeasible );
         SCIPdebug( SCIP_CALL( SCIPprintRow(scip, rowub, NULL) ) );
         ++tcliquedata->nboundcuts;
         ++tcliquedata->ncuts;
         *success = TRUE;
      }
      SCIP_CALL( SCIPreleaseRow(scip, &rowub) );
   }

   return SCIP_OKAY;
}


/** Generate bound constraint
 *
 *  We generate the row corresponding to the following simple valid inequalities:
 *  \f[
 *         \frac{x_1}{u_1} + \ldots + \frac{x_n}{u_n} \leq 1\qquad\mbox{and}\qquad
 *         \frac{x_1}{\ell_1} + \ldots + \frac{x_n}{\ell_1} \leq 1,
 *  \f]
 *  where \f$\ell_1, \ldots, \ell_n\f$ and \f$u_1, \ldots, u_n\f$ are the nonzero and finite lower and upper bounds of
 *  the variables \f$x_1, \ldots, x_n\f$. If an upper bound < 0 or a lower bound > 0, the constraint itself is
 *  redundant, so the cut is not applied (lower bounds > 0 and upper bounds < 0 are usually detected in presolving or
 *  propagation). Infinite bounds and zero are skipped. Thus \f$\ell_1, \ldots, \ell_n\f$ are all negative, which
 *  results in the \f$\leq\f$ inequality. In case of the presence of variable upper bounds, the bound inequality can
 *  be further strengthened.
 *
 *  Note that in fact, any mixture of nonzero finite lower and upper bounds would lead to a valid inequality as
 *  above. However, usually either the lower or upper bound is nonzero. Thus, the above inequalities are the most
 *  interesting.
 */
static
SCIP_RETCODE SCIPgenerateBoundInequalityFromSOS1Nodes(
   SCIP*                 scip,               /**< SCIP pointer */
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler */
   SCIP_DIGRAPH*         conflictgraph,      /**< conflict graph */
   int*                  nodes,              /**< conflict graph nodes for bound constraint */
   int                   nnodes,             /**< number of conflict graph nodes for bound constraint */
   SCIP_Real             rhs,                /**< right hand side of bound constraint */
   SCIP_Bool             local,              /**< in any case produce a local cut (even if local bounds of variables are valid globally) */
   SCIP_Bool             global,             /**< in any case produce a global cut */
   SCIP_Bool             strengthen,         /**< whether trying to strengthen bound constraint */
   SCIP_Bool             removable,          /**< should the inequality be removed from the LP due to aging or cleanup? */
   const char *          nameext,            /**< part of name of bound constraints */
   SCIP_ROW**            rowlb,              /**< output: row for lower bounds (or NULL if not needed) */
   SCIP_ROW**            rowub               /**< output: row for upper bounds (or NULL if not needed) */
   )
{
   char name[SCIP_MAXSTRLEN];
   SCIP_VAR* lbboundvar = NULL;
   SCIP_VAR* ubboundvar = NULL;
   SCIP_Bool locallbs;
   SCIP_Bool localubs;
   SCIP_VAR** vars;
   SCIP_Real* vals;

   assert( scip != NULL );
   assert( conshdlr != NULL );
   assert( conflictgraph != NULL );
   assert( ! local || ! global );
   assert( nodes != NULL );

   /* allocate buffer array */
   SCIP_CALL( SCIPallocBufferArray(scip, &vars, nnodes+1) );
   SCIP_CALL( SCIPallocBufferArray(scip, &vals, nnodes+1) );

   /* take care of upper bounds */
   if ( rowub != NULL )
   {
      SCIP_Bool useboundvar;
      int cnt;
      int j;

      /* loop through all variables. We check whether all bound variables (if existent) are equal; if this is the
       * case then the bound constraint can be strengthened */
      cnt = 0;
      localubs = local;
      useboundvar = strengthen;
      for (j = 0; j < nnodes; ++j)
      {
         SCIP_NODEDATA* nodedata;
         SCIP_VAR* var;
         SCIP_Real val;

         nodedata = (SCIP_NODEDATA*)SCIPdigraphGetNodeData(conflictgraph, nodes[j]);
         assert( nodedata != NULL );
         var = nodedata->var;
         assert( var != NULL );

         /* if variable is not involved in a variable bound constraint */
         if ( ! useboundvar || nodedata->ubboundvar == NULL )
         {
            useboundvar = FALSE;
            if ( localubs )
            {
               assert( ! global );
               val = SCIPvarGetUbLocal(var);
            }
            else
            {
               val = SCIPvarGetUbGlobal(var);
               if ( ! global && ! SCIPisFeasEQ(scip, val, SCIPvarGetUbLocal(var)) )
               {
                  localubs = TRUE;

                  /* restart 'for'-loop, since we need the local bounds of the variables */
                  j = -1;
                  cnt = 0;
                  continue;
               }
            }
         }
         else
         {
            /* in this case the cut is always valid globally */

            /* if we have a bound variable for the first time */
            if ( ubboundvar == NULL )
            {
               ubboundvar = nodedata->ubboundvar;
               val = nodedata->ubboundcoef;
            }
            /* else if the bound variable equals the stored bound variable */
            else if ( SCIPvarCompare(ubboundvar, nodedata->ubboundvar) == 0 )
            {
               val = nodedata->ubboundcoef;
            }
            else /* else use bounds on the variables */
            {
               useboundvar = FALSE;

               /* restart 'for'-loop */
               j = -1;
               cnt = 0;
               continue;
            }
         }

         /* should not apply the cut if a variable is fixed to be negative -> constraint is redundant */
         if ( SCIPisNegative(scip, val) )
            break;

         /* store variable if relevant for bound inequality */
         if ( ! SCIPisInfinity(scip, val) && ! SCIPisZero(scip, val) )
         {
            vars[cnt] = var;
            vals[cnt++] = 1.0/val;
         }
      }

      /* if cut is meaningful */
      if ( j == nnodes && cnt >= 2 )
      {
         if ( useboundvar )
         {
            /* add bound variable to array */
            vars[cnt] = ubboundvar;
            vals[cnt++] = -rhs;
            assert(ubboundvar != NULL );

            /* create upper bound inequality if at least two of the bounds are finite and nonzero */
            (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "sosub#%s", nameext);
            SCIP_CALL( SCIPcreateEmptyRowCons(scip, rowub, conshdlr, name, -SCIPinfinity(scip), 0.0, localubs, FALSE, removable) );
            SCIP_CALL( SCIPaddVarsToRow(scip, *rowub, cnt, vars, vals) );
            SCIPdebug( SCIP_CALL( SCIPprintRow(scip, *rowub, NULL) ) );
         }
         else
         {
            /* create upper bound inequality if at least two of the bounds are finite and nonzero */
            (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "sosub#%s", nameext);
            SCIP_CALL( SCIPcreateEmptyRowCons(scip, rowub, conshdlr, name, -SCIPinfinity(scip), rhs, localubs, FALSE, removable) );
            SCIP_CALL( SCIPaddVarsToRow(scip, *rowub, cnt, vars, vals) );
            SCIPdebug( SCIP_CALL( SCIPprintRow(scip, *rowub, NULL) ) );
         }
      }
   }


   /* take care of lower bounds */
   if ( rowlb != NULL )
   {
      SCIP_Bool useboundvar;
      int cnt;
      int j;

      /* loop through all variables. We check whether all bound variables (if existent) are equal; if this is the
       * case then the bound constraint can be strengthened */
      cnt = 0;
      locallbs = local;
      useboundvar = strengthen;
      for (j = 0; j < nnodes; ++j)
      {
         SCIP_NODEDATA* nodedata;
         SCIP_VAR* var;
         SCIP_Real val;

         nodedata = (SCIP_NODEDATA*)SCIPdigraphGetNodeData(conflictgraph, nodes[j]);
         assert( nodedata != NULL );
         var = nodedata->var;
         assert( var != NULL );

         /* if variable is not involved in a variable bound constraint */
         if ( ! useboundvar || nodedata->lbboundvar == NULL )
         {
            useboundvar = FALSE;
            if ( locallbs )
            {
               assert( ! global );
               val = SCIPvarGetLbLocal(var);
            }
            else
            {
               val = SCIPvarGetLbGlobal(var);
               if ( ! global && ! SCIPisFeasEQ(scip, val, SCIPvarGetLbLocal(var)) )
               {
                  locallbs = TRUE;

                  /* restart 'for'-loop, since we need the local bounds of the variables */
                  j = -1;
                  cnt = 0;
                  continue;
               }
            }
         }
         else
         {
            /* in this case the cut is always valid globally */

            /* if we have a bound variable for the first time */
            if ( lbboundvar == NULL )
            {
               lbboundvar = nodedata->lbboundvar;
               val = nodedata->lbboundcoef;
            }
            /* else if the bound variable equals the stored bound variable */
            else if ( SCIPvarCompare(lbboundvar, nodedata->lbboundvar) == 0 )
            {
               val = nodedata->lbboundcoef;
            }
            else /* else use bounds on the variables */
            {
               useboundvar = FALSE;

               /* restart 'for'-loop */
               j = -1;
               cnt = 0;
               continue;
            }
         }

         /* should not apply the cut if a variable is fixed to be positive -> constraint is redundant */
         if ( SCIPisPositive(scip, val) )
            break;

         /* store variable if relevant for bound inequality */
         if ( ! SCIPisInfinity(scip, val) && ! SCIPisZero(scip, val) )
         {
            vars[cnt] = var;
            vals[cnt++] = 1.0/val;
         }
      }

      /* if cut is meaningful */
      if ( j == nnodes && cnt >= 2 )
      {
         if ( useboundvar )
         {
            /* add bound variable to array */
            vars[cnt] = lbboundvar;
            vals[cnt++] = -rhs;
            assert(lbboundvar != NULL );

            /* create upper bound inequality if at least two of the bounds are finite and nonzero */
            (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "soslb#%s", nameext);
            SCIP_CALL( SCIPcreateEmptyRowCons(scip, rowlb, conshdlr, name, -SCIPinfinity(scip), 0.0, locallbs, FALSE, TRUE) );
            SCIP_CALL( SCIPaddVarsToRow(scip, *rowlb, cnt, vars, vals) );
            SCIPdebug( SCIP_CALL( SCIPprintRow(scip, *rowlb, NULL) ) );
         }
         else
         {
            /* create upper bound inequality if at least two of the bounds are finite and nonzero */
            (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "soslb#%s", nameext);
            SCIP_CALL( SCIPcreateEmptyRowCons(scip, rowlb, conshdlr, name, -SCIPinfinity(scip), rhs, locallbs, FALSE, TRUE) );
            SCIP_CALL( SCIPaddVarsToRow(scip, *rowlb, cnt, vars, vals) );
            SCIPdebug( SCIP_CALL( SCIPprintRow(scip, *rowlb, NULL) ) );
         }
      }
   }

   /* free buffer array */
   SCIPfreeBufferArray(scip, &vals);
   SCIPfreeBufferArray(scip, &vars);

   return SCIP_OKAY;
}


/** generates bound cuts using a clique found by algorithm for maximum weight clique
 *  and decides whether to stop generating cliques with the algorithm for maximum weight clique
 */
static
TCLIQUE_NEWSOL(tcliqueNewsolClique)
{
   TCLIQUE_WEIGHT minweightinc;

   assert( acceptsol != NULL );
   assert( stopsolving != NULL );
   assert( tcliquedata != NULL );

   /* we don't accept the solution as new incumbent, because we want to find many violated clique inequalities */
   *acceptsol = FALSE;
   *stopsolving = FALSE;

   /* slightly increase the minimal weight for additional cliques */
   minweightinc = (cliqueweight - *minweight)/10;
   minweightinc = MAX(minweightinc, 1);
   *minweight += minweightinc;

   /* adds cut if weight of the clique is greater than 1 */
   if( cliqueweight > tcliquedata->scaleval )
   {
      SCIP* scip;
      SCIP_SOL* sol;
      SCIP_Real unscaledweight;
      SCIP_Real solval;
      SCIP_Real bound;
      SCIP_VAR* var;
      int node;
      int i;

      scip = tcliquedata->scip;
      sol = tcliquedata->sol;
      assert( scip != NULL );

      /* calculate the weight of the clique in unscaled fractional variable space */
      unscaledweight = 0.0;
      for( i = 0; i < ncliquenodes; i++ )
      {
         node = cliquenodes[i];
         var = nodeGetVarSOS1(tcliquedata->conflictgraph, node);
         solval = SCIPgetSolVal(scip, sol, var);

         if ( SCIPisFeasPositive(scip, solval) )
         {
            if ( tcliquedata->strthenboundcuts )
               bound = REALABS( SCIPnodeGetSolvalVarboundUbSOS1(scip, tcliquedata->conflictgraph, sol, node) );
            else
               bound = REALABS( SCIPvarGetUbLocal(var) );
         }
         else if ( SCIPisFeasNegative(scip, solval) )
         {
            if ( tcliquedata->strthenboundcuts )
               bound = REALABS( SCIPnodeGetSolvalVarboundLbSOS1(scip, tcliquedata->conflictgraph, sol, node) );
            else
               bound = REALABS( SCIPvarGetLbLocal(var) );
         }
         else
            bound = 0.0;

         solval = REALABS( solval );

         if ( ! SCIPisFeasZero(scip, bound) && ! SCIPisInfinity(scip, bound) )
            unscaledweight += REALABS( solval/bound );
      }

      if( SCIPisEfficacious(scip, unscaledweight - 1.0) )
      {
         char nameext[SCIP_MAXSTRLEN];
         SCIP_ROW* rowlb = NULL;
         SCIP_ROW* rowub = NULL;
         SCIP_Bool success;

         /* generate bound inequalities for lower and upper bound case
          * NOTE: tests have shown that non-removable rows give the best results */
         (void) SCIPsnprintf(nameext, SCIP_MAXSTRLEN, "%d", tcliquedata->nboundcuts);
         if( SCIPgenerateBoundInequalityFromSOS1Nodes(scip, tcliquedata->conshdlr, tcliquedata->conflictgraph,
               cliquenodes, ncliquenodes, 1.0, FALSE, FALSE, tcliquedata->strthenboundcuts, FALSE, nameext, &rowlb, &rowub) != SCIP_OKAY )
         {
            SCIPerrorMessage("unexpected error in bound cut creation.\n");
            SCIPABORT();
         }

         /* add bound cut(s) to separation storage if existent */
         if ( addBoundCutSepa(scip, tcliquedata, rowlb, rowub, &success) != SCIP_OKAY )
         {
            SCIPerrorMessage("unexpected error in bound cut creation.\n");
            SCIPABORT();
         }

         /* if at least one cut has been added */
         if ( success )
         {
            SCIPdebugMessage(" -> found bound cut corresponding to clique (act=%g)\n", unscaledweight);

            /* if we found more than half the cuts we are allowed to generate, we accept the clique as new incumbent,
             * such that only more violated cuts are generated afterwards
             */
            if( tcliquedata->maxboundcuts >= 0 )
            {
               if( tcliquedata->ncuts > tcliquedata->maxboundcuts/2 )
                  *acceptsol = TRUE;
               if( tcliquedata->ncuts >= tcliquedata->maxboundcuts )
                  *stopsolving = TRUE;
            }
         }
         else
            *stopsolving = TRUE;
      }
   }
}


/** separate bound inequalities from conflict graph */
static
SCIP_RETCODE sepaBoundInequalitiesFromGraph(
   SCIP*                 scip,               /**< SCIP pointer */
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler */
   SCIP_CONSHDLRDATA*    conshdlrdata,       /**< constraint handler data */
   SCIP_SOL*             sol,                /**< LP solution to be separated (or NULL) */
   int                   maxboundcuts,       /**< maximal number of bound cuts separated per separation round (-1: no limit) */
   int*                  ngen,               /**< pointer to store number of cuts generated */
   SCIP_RESULT*          result              /**< pointer to store result of separation */
   )
{
   SCIP_DIGRAPH* conflictgraph;
   TCLIQUE_DATA* tcliquedata;
   TCLIQUE_WEIGHT cliqueweight;
   TCLIQUE_STATUS tcliquestatus;
   int nsos1vars;

   SCIP_Real scaleval = 1000.0;                  /* factor for scaling weights */
   int maxtreenodes = 10000;                     /* maximal number of nodes of b&b tree */
   int maxzeroextensions = 1000;                 /* maximal number of zero-valued variables extending the clique (-1: no limit) */
   int backtrackfreq = 1000;                     /* frequency for premature backtracking up to tree level 1 (0: no backtracking) */
   int ntreenodes;
   int* cliquenodes;
   int ncliquenodes;

   assert( scip != NULL );
   assert( conshdlr != NULL );
   assert( conshdlrdata != NULL );
   assert( ngen != NULL );
   assert( result != NULL );

   /* get conflict graph */
   conflictgraph = SCIPgetConflictgraphSOS1(conshdlr);
   assert( conflictgraph != NULL );

   /* get number of SOS1 variables */
   nsos1vars = SCIPgetNSOS1Vars(conshdlr);

   /* initialize tclique graph if not done already */
   if ( conshdlrdata->tcliquegraph == NULL )
   {
      SCIP_CALL( initTCliquegraph(scip, conshdlr, conshdlrdata, conflictgraph, nsos1vars, sol, scaleval) );
   }
   tcliquedata = conshdlrdata->tcliquedata;
   tcliquedata->maxboundcuts = maxboundcuts;
   tcliquedata->ncuts = 0;

   /* update the weights of the tclique graph */
   SCIP_CALL( updateWeightsTCliquegraph(scip, conshdlrdata, tcliquedata, conflictgraph, sol, nsos1vars) );

   /* allocate buffer array */
   SCIP_CALL( SCIPallocBufferArray(scip, &cliquenodes, nsos1vars) );

   /* start algorithm to find maximum weight cliques and use them to generate bound cuts */
   tcliqueMaxClique(tcliqueGetNNodes, tcliqueGetWeights, tcliqueIsEdge, tcliqueSelectAdjnodes,
      conshdlrdata->tcliquegraph, tcliqueNewsolClique, tcliquedata,
      cliquenodes, &ncliquenodes, &cliqueweight, (int)scaleval-1, (int)scaleval+1,
      maxtreenodes, backtrackfreq, maxzeroextensions, -1, &ntreenodes, &tcliquestatus);

   /* free buffer array */
   SCIPfreeBufferArray(scip, &cliquenodes);

   /* get number of cuts of current separation round */
   *ngen = tcliquedata->ncuts;

   /* update number of bound cuts in separator data */
   conshdlrdata->nboundcuts = tcliquedata->nboundcuts;

   /* evaluate the result of the separation */
   if ( *ngen > 0 )
      *result = SCIP_SEPARATED;
   else
      *result = SCIP_DIDNOTFIND;

   return SCIP_OKAY;
}


/** Generate a bound constraint from the variables of an SOS1 constraint (see SCIPgenerateBoundInequalityFromSOS1Nodes() for more information) */
static
SCIP_RETCODE SCIPgenerateBoundInequalityFromSOS1Cons(
   SCIP*                 scip,               /**< SCIP pointer */
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler */
   SCIP_CONS*            cons,               /**< SOS1 constraint */
   SCIP_Bool             local,              /**< in any case produce a local cut (even if local bounds of variables are valid globally) */
   SCIP_Bool             global,             /**< in any case produce a global cut */
   SCIP_Bool             strengthen,         /**< whether trying to strengthen bound constraint */
   SCIP_Bool             removable,          /**< should the inequality be removed from the LP due to aging or cleanup? */
   SCIP_ROW**            rowlb,              /**< output: row for lower bounds (or NULL if not needed) */
   SCIP_ROW**            rowub               /**< output: row for upper bounds (or NULL if not needed) */
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONSDATA* consdata;
   int* nodes;
   int nvars;
   int j;

   assert( scip != NULL );
   assert( conshdlr != NULL );
   assert( cons != NULL );

   /* get constraint data */
   consdata = SCIPconsGetData(cons);
   assert( consdata != NULL );
   assert( consdata->vars != NULL );
   nvars = consdata->nvars;

   /* get constraint handler data */
   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert( conshdlrdata != NULL );
   assert( strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0 );

   /* allocate buffer array */
   SCIP_CALL( SCIPallocBufferArray(scip, &nodes, nvars) );

   /* get nodes in the conflict graph */
   for (j = 0; j < nvars; ++j)
   {
      nodes[j] = varGetNodeSOS1(conshdlr, consdata->vars[j]);
      assert( nodes[j] >= 0 );
   }

   /* generate bound constraint from conflict graph nodes */
   SCIP_CALL( SCIPgenerateBoundInequalityFromSOS1Nodes(scip, conshdlr, conshdlrdata->conflictgraph, nodes, nvars, 1.0, local, global, strengthen, removable, SCIPconsGetName(cons), rowlb, rowub) );

   /* free buffer array */
   SCIPfreeBufferArray(scip, &nodes);

   return SCIP_OKAY;
}


/** initialize or separate bound inequalities from SOS1 constraints */
static
SCIP_RETCODE initsepaBoundInequalityFromSOS1Cons(
   SCIP*                 scip,               /**< SCIP pointer */
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler */
   SCIP_CONSHDLRDATA*    conshdlrdata,       /**< constraint handler data */
   SCIP_CONS**           conss,              /**< SOS1 constraints */
   int                   nconss,             /**< number of SOS1 constraints */
   SCIP_SOL*             sol,                /**< LP solution to be separated (or NULL) */
   SCIP_Bool             solvedinitlp,       /**< TRUE if initial LP relaxation at a node is solved */
   int                   maxboundcuts,       /**< maximal number of bound cuts separated per separation round (-1: no limit) */
   int*                  ngen,               /**< pointer to store number of cuts generated (or NULL) */
   SCIP_RESULT*          result              /**< pointer to store result of separation (or NULL) */
   )
{
   SCIP_Bool cutoff = FALSE;
   int c;

   assert( scip != NULL );
   assert( conshdlrdata != NULL );
   assert( conss != NULL );

   if ( result != NULL )
      *result = SCIP_DIDNOTFIND;

   for (c = 0; c < nconss; ++c)
   {
      SCIP_CONSDATA* consdata;
      SCIP_ROW* row;

      assert( conss != NULL );
      assert( conss[c] != NULL );
      consdata = SCIPconsGetData(conss[c]);
      assert( consdata != NULL );

      if ( solvedinitlp )
         SCIPdebugMessage("Separating inequalities for SOS1 constraint <%s>.\n", SCIPconsGetName(conss[c]) );
      else
         SCIPdebugMessage("Checking for initial rows for SOS1 constraint <%s>.\n", SCIPconsGetName(conss[c]) );

      /* possibly generate rows if not yet done */
      if ( consdata->rowub == NULL || consdata->rowlb == NULL )
      {
         SCIP_ROW* rowlb = NULL;
         SCIP_ROW* rowub = NULL;

         SCIP_CALL( SCIPgenerateBoundInequalityFromSOS1Cons(scip, conshdlr, conss[c], FALSE, TRUE, TRUE, FALSE, &rowlb, &rowub) );

         /* if row(s) should be globally stored in constraint data */
         if ( rowlb != NULL )
         {
            consdata->rowlb = rowlb;
         }
         if ( rowub != NULL )
         {
            consdata->rowub = rowub;
         }
      }

      /* put corresponding rows into LP */
      row = consdata->rowub;
      if ( row != NULL && ! SCIProwIsInLP(row) && ( solvedinitlp || SCIPisCutEfficacious(scip, sol, row) ) )
      {
         assert( SCIPisInfinity(scip, -SCIProwGetLhs(row)) && ( SCIPisEQ(scip, SCIProwGetRhs(row), 1.0) || SCIPisEQ(scip, SCIProwGetRhs(row), 0.0) ) );

         SCIP_CALL( SCIPaddCut(scip, NULL, row, FALSE, &cutoff) );
         if ( cutoff && result != NULL )
            break;
         assert( ! cutoff );
         SCIPdebug( SCIP_CALL( SCIPprintRow(scip, row, NULL) ) );

         if ( solvedinitlp )
         {
            assert( ngen != NULL );
            SCIP_CALL( SCIPresetConsAge(scip, conss[c]) );
            ++(*ngen);
         }
      }
      row = consdata->rowlb;
      if ( row != NULL && ! SCIProwIsInLP(row) && ( solvedinitlp || SCIPisCutEfficacious(scip, sol, row) ) )
      {
         assert( SCIPisInfinity(scip, -SCIProwGetLhs(row)) && ( SCIPisEQ(scip, SCIProwGetRhs(row), 1.0) || SCIPisEQ(scip, SCIProwGetRhs(row), 0.0) ) );

         SCIP_CALL( SCIPaddCut(scip, NULL, row, FALSE, &cutoff) );
         if ( cutoff && result != NULL )
            break;
         assert( ! cutoff );
         SCIPdebug( SCIP_CALL( SCIPprintRow(scip, row, NULL) ) );

         if ( solvedinitlp )
         {
            assert( ngen != NULL );
            SCIP_CALL( SCIPresetConsAge(scip, conss[c]) );
            ++(*ngen);
         }
      }

      if ( ngen != NULL && maxboundcuts >= 0 && *ngen >= maxboundcuts )
         break;
   }

   if ( cutoff )
      *result = SCIP_CUTOFF;
   else if ( ngen != NULL && *ngen > 0 )
      *result = SCIP_SEPARATED;

   return SCIP_OKAY;
}


/* --------------------initialization/deinitialization ------------------------*/

/** check whether var1 is a bound variable of var0; i.e., var0 >= c * var1 or var0 <= d * var1.
 *  If true, then add this information to the node data of the conflict graph.
 */
static
SCIP_RETCODE detectVarboundSOS1(
   SCIP*                 scip,               /**< SCIP pointer */
   SCIP_CONSHDLR*        conshdlr,           /**< SOS1 constraint handler */
   SCIP_CONSHDLRDATA*    conshdlrdata,       /**< SOS1 constraint handler data */
   SCIP_VAR*             var0,               /**< first variable */
   SCIP_VAR*             var1,               /**< second variable */
   SCIP_Real             val0,               /**< first coefficient */
   SCIP_Real             val1                /**< second coefficient */
   )
{
   int node0;

   assert( scip != NULL );
   assert( conshdlr != NULL );
   assert( conshdlrdata != NULL );
   assert( var0 != NULL && var1 != NULL );

   /* get nodes of variable in the conflict graph (node = -1 if no SOS1 variable) */
   node0 = varGetNodeSOS1(conshdlr, var0);

   /* if var0 is an SOS1 variable */
   if ( node0 >= 0 )
   {
      SCIP_Real val;

      assert( ! SCIPisFeasZero(scip, val0) );
      val = -val1/val0;

      /* check variable bound relation of variables */

      /* handle lower bound case */
      if ( SCIPisFeasNegative(scip, val0) && SCIPisFeasNegative(scip, val) )
      {
         SCIP_NODEDATA* nodedata;

         /* get node data of the conflict graph */
         nodedata = (SCIP_NODEDATA*)SCIPdigraphGetNodeData(conshdlrdata->conflictgraph, node0);

         /* @todo: maybe save multiple variable bounds for each SOS1 variable */
         if ( nodedata->lbboundvar == NULL )
         {
            /* add variable bound information to node data */
            nodedata->lbboundvar = var1;
            nodedata->lbboundcoef = val;

            SCIPdebugMessage("detected variable bound constraint %s >= %f %s.\n", SCIPvarGetName(var0), val, SCIPvarGetName(var1));
         }
      }
      /* handle upper bound case */
      else if ( SCIPisFeasPositive(scip, val0) && SCIPisFeasPositive(scip, val) )
      {
         SCIP_NODEDATA* nodedata;
         assert( SCIPisFeasPositive(scip, val0) );

         /* get node data of the conflict graph */
         nodedata = (SCIP_NODEDATA*)SCIPdigraphGetNodeData(conshdlrdata->conflictgraph, node0);

         if ( nodedata->ubboundvar == NULL )
         {
            /* add variable bound information to node data */
            nodedata->ubboundvar = var1;
            nodedata->ubboundcoef = val;

            SCIPdebugMessage("detected variable bound constraint %s <= %f %s.\n", SCIPvarGetName(var0), val, SCIPvarGetName(var1));
         }
      }
   }

   return SCIP_OKAY;
}


/* pass connected component @p C of the conflict graph and check whether all the variables correspond to a unique variable upper bound variable @p z,
 *  i.e., \f$x_i \leq u_i z\f$ for every \f$i\in C\f$.
 *
 *  Note: if the upper bound variable is not unique, then bound inequalities usually cannot be strengthened.
 */
static
SCIP_RETCODE passConComponentVarbound(
   SCIP*                 scip,               /**< SCIP pointer */
   SCIP_DIGRAPH*         conflictgraph,      /**< conflict graph */
   int                   node,               /**< current node of connected component */
   SCIP_VAR*             boundvar,           /**< bound variable of connected component */
   SCIP_Bool             checklb,            /**< whether to check lower bound variable (else upper bound variable) */
   SCIP_Bool*            processed,          /**< states for each variable whether it has been processed */
   int*                  concomp,            /**< current connected component */
   int*                  nconcomp,           /**< pointer to store number of elements of connected component */
   SCIP_Bool*            unique              /**< pointer to store whether bound variable is unique */
   )
{
   int* succ;
   int nsucc;
   int s;

   assert( scip != NULL );
   assert( conflictgraph != NULL );
   assert( processed != NULL );
   assert( concomp != NULL );
   assert( nconcomp != NULL );
   assert( unique != NULL );

   processed[node] = TRUE;
   concomp[(*nconcomp)++] = node;

   /* if bound variable of connected component without new node is unique */
   if ( unique )
   {
      SCIP_NODEDATA* nodedata;
      SCIP_VAR* comparevar;
      nodedata = (SCIP_NODEDATA*)SCIPdigraphGetNodeData(conflictgraph, node);
      assert( nodedata != NULL );

      if ( checklb )
         comparevar = nodedata->lbboundvar;
      else
         comparevar = nodedata->ubboundvar;

      /* check whether bound variable is unique for connected component without new node */
      if ( boundvar == NULL )
      {
         if ( comparevar != NULL )
            unique = FALSE;
      }
      else
      {
         if ( comparevar == NULL )
            unique = FALSE;
         else if ( SCIPvarCompare(boundvar, comparevar) != 0 )
            unique = FALSE;
      }
   }

   /* pass through successor variables */
   nsucc = SCIPdigraphGetNSuccessors(conflictgraph, node);
   succ = SCIPdigraphGetSuccessors(conflictgraph, node);
   for (s = 0; s < nsucc; ++s)
   {
      if ( ! processed[succ[s]] )
         SCIP_CALL( passConComponentVarbound(scip, conflictgraph, succ[s], boundvar, checklb, processed, concomp, nconcomp, unique) );
   }

   return SCIP_OKAY;
}


/** for each connected component @p C of the conflict graph check whether all the variables correspond to a unique variable upper bound variable @p z
 *  (e.g., for the upper bound case this means that \f$x_i \leq u_i z\f$ for every \f$i\in C\f$).
 *
 *  Note: if the bound variable is not unique, then bound inequalities usually cannot be strengthened.
 */
static
SCIP_RETCODE checkConComponentsVarbound(
   SCIP*                 scip,               /**< SCIP pointer */
   SCIP_DIGRAPH*         conflictgraph,      /**< conflict graph */
   int                   nsos1vars,          /**< number of SOS1 variables */
   SCIP_Bool             checklb             /**< whether to check lower bound variable (else check upper bound variable) */
   )
{
   SCIP_Bool* processed;  /* states for each variable whether it has been processed */
   int* concomp;          /* current connected component */
   int nconcomp;
   int j;

   assert( scip != NULL );
   assert( conflictgraph != NULL );

   /* allocate buffer arrays and initialize 'processed' array */
   SCIP_CALL( SCIPallocBufferArray(scip, &processed, nsos1vars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &concomp, nsos1vars) );
   for (j = 0; j < nsos1vars; ++j)
      processed[j] = FALSE;

   /* run through all SOS1 variables */
   for (j = 0; j < nsos1vars; ++j)
   {
      /* if variable belongs to a connected component that has not been processed so far */
      if ( ! processed[j] )
      {
         SCIP_NODEDATA* nodedata;
         SCIP_VAR* boundvar;
         SCIP_Bool unique;
         int* succ;
         int nsucc;
         int s;

         nodedata = (SCIP_NODEDATA*)SCIPdigraphGetNodeData(conflictgraph, j);
         assert( nodedata != NULL );

         if ( checklb )
            boundvar = nodedata->lbboundvar;
         else
            boundvar = nodedata->ubboundvar;
         unique = TRUE;

         processed[j] = TRUE;
         concomp[0] = j;
         nconcomp = 1;

         /* pass through successor variables */
         nsucc = SCIPdigraphGetNSuccessors(conflictgraph, j);
         succ = SCIPdigraphGetSuccessors(conflictgraph, j);
         for (s = 0; s < nsucc; ++s)
         {
            if ( ! processed[succ[s]] )
               SCIP_CALL( passConComponentVarbound(scip, conflictgraph, succ[s], boundvar, checklb, processed, concomp, &nconcomp, &unique) );
         }

         /* if the connected component has a unique bound variable */
         if ( unique && boundvar != NULL )
         {
            for (s = 0; s < nconcomp; ++s)
            {
               nodedata = (SCIP_NODEDATA*)SCIPdigraphGetNodeData(conflictgraph, concomp[s]);
               assert( processed[concomp[s]] == TRUE );
               assert( nodedata != NULL );

               if ( checklb )
                  nodedata->lbboundcomp = TRUE;
               else
                  nodedata->ubboundcomp = TRUE;
            }
            SCIPdebugMessage("Found a connected component of size <%i> with unique bound variable.\n", nconcomp);
         }
      }
   }

   /* free buffer arrays */
   SCIPfreeBufferArray(scip, &concomp);
   SCIPfreeBufferArray(scip, &processed);

   return SCIP_OKAY;
}


/** check all linear constraints for variable bound constraints of the form c*z <= x <= d*z, where @p x is some SOS1
 *  variable and @p z some arbitrary variable (not necessarily binary)
 */
static
SCIP_RETCODE checkLinearConssVarboundSOS1(
   SCIP*                 scip,               /**< SCIP pointer */
   SCIP_CONSHDLR*        conshdlr,           /**< SOS1 constraint handler */
   SCIP_CONSHDLRDATA*    conshdlrdata,       /**< SOS1 constraint handler data */
   SCIP_CONS**           linconss,           /**< linear constraints */
   int                   nlinconss           /**< number of linear constraints */
   )
{
   int c;

   /* loop through linear constraints */
   for (c = 0; c < nlinconss; ++c)
   {
      SCIP_CONS* lincons;
      int nvars;

      lincons = linconss[c];

      /* variable bound constraints only contain two variables */
      nvars = SCIPgetNVarsLinear(scip, lincons);
      if ( nvars == 2 )
      {
         SCIP_VAR** vars;
         SCIP_Real* vals;
         SCIP_VAR* var0;
         SCIP_VAR* var1;
         SCIP_Real lhs;
         SCIP_Real rhs;

         /* get constraint data */
         vars = SCIPgetVarsLinear(scip, lincons);
         vals = SCIPgetValsLinear(scip, lincons);
         lhs = SCIPgetLhsLinear(scip, lincons);
         rhs = SCIPgetRhsLinear(scip, lincons);

         var0 = vars[0];
         var1 = vars[1];
         assert( var0 != NULL && var1 != NULL );

         /* at least one variable should be an SOS1 variable */
         if ( varIsSOS1(conshdlr, var0) || varIsSOS1(conshdlr, var1) )
         {
            SCIP_Real val0;
            SCIP_Real val1;

            /* check whether right hand side or left hand side of constraint is zero */
            if ( SCIPisFeasZero(scip, lhs) )
            {
               val0 = -vals[0];
               val1 = -vals[1];

               /* check whether the two variables are in a variable bound relation */
               SCIP_CALL( detectVarboundSOS1(scip, conshdlr, conshdlrdata, var0, var1, val0, val1) );
               SCIP_CALL( detectVarboundSOS1(scip, conshdlr, conshdlrdata, var1, var0, val1, val0) );
            }
            else if( SCIPisFeasZero(scip, rhs) )
            {
               val0 = vals[0];
               val1 = vals[1];

               /* check whether the two variables are in a variable bound relation */
               SCIP_CALL( detectVarboundSOS1(scip, conshdlr, conshdlrdata, var0, var1, val0, val1) );
               SCIP_CALL( detectVarboundSOS1(scip, conshdlr, conshdlrdata, var1, var0, val1, val0) );
            }
         }
      }
   }

   return SCIP_OKAY;
}


/** set node data of conflict graph nodes */
static
SCIP_RETCODE setNodeDataSOS1(
   SCIP*                 scip,               /**< SCIP pointer */
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler */
   SCIP_CONSHDLRDATA*    conshdlrdata,       /**< SOS1 constraint handler data */
   int                   nsos1conss,         /**< number of SOS1 constraints */
   int                   nsos1vars           /**< number of SOS1 variables */
   )
{
   SCIP_CONSHDLR* linconshdlr;
   SCIP_CONS** linconss;
   int nlinconss;

   /* if no SOS1 variables exist -> exit */
   if ( nsos1vars == 0 )
      return SCIP_OKAY;

   /* get constraint handler data of linear constraints */
   linconshdlr = SCIPfindConshdlr(scip, "linear");
   if ( linconshdlr == NULL )
      return SCIP_OKAY;

   /* get linear constraints and number of linear constraints */
   nlinconss = SCIPconshdlrGetNConss(linconshdlr);
   linconss = SCIPconshdlrGetConss(linconshdlr);

   /* check linear constraints for variable bound constraints */
   SCIP_CALL( checkLinearConssVarboundSOS1(scip, conshdlr, conshdlrdata, linconss, nlinconss) );

   /* for each connected component of the conflict graph check whether all the variables correspond to a unique variable
    * upper bound variable */
   SCIP_CALL( checkConComponentsVarbound(scip, conshdlrdata->conflictgraph, conshdlrdata->nsos1vars, TRUE) );
   SCIP_CALL( checkConComponentsVarbound(scip, conshdlrdata->conflictgraph, conshdlrdata->nsos1vars, FALSE) );

   return SCIP_OKAY;
}


/* initialize conflictgraph and create hashmap for SOS1 variables */
static
SCIP_RETCODE initConflictgraph(
   SCIP*                 scip,               /**< SCIP pointer */
   SCIP_CONSHDLRDATA*    conshdlrdata,       /**< constraint handler data */
   SCIP_CONS**           conss,              /**< SOS1 constraints */
   int                   nconss              /**< number of SOS1 constraints */
   )
{
   SCIP_Bool* nodecreated; /* nodecreated[i] = TRUE if a node in the conflictgraph is already created for index i
                            * (with i index of the original variables) */
   int* nodeorig;          /* nodeorig[i] = node of original variable x_i in the conflictgraph */
   int ntotalvars;
   int cntsos;
   int i;
   int j;
   int c;

   assert( conshdlrdata != NULL );
   assert( nconss == 0 || conss != NULL );

   /* get the number of original problem variables */
   ntotalvars = SCIPgetNTotalVars(scip);

   /* initialize vector 'nodecreated' */
   SCIP_CALL( SCIPallocBufferArray(scip, &nodeorig, ntotalvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &nodecreated, ntotalvars) );
   for (i = 0; i < ntotalvars; ++i)
      nodecreated[i] = FALSE;

   /* compute number of SOS1 variables */
   cntsos = 0;
   for (c = 0; c < nconss; ++c)
   {
      SCIP_CONSDATA* consdata;
      SCIP_VAR** vars;
      int nvars;

      assert( conss[c] != NULL );

      /* get constraint data field of the constraint */
      consdata = SCIPconsGetData(conss[c]);
      assert( consdata != NULL );

      /* get variables and number of variables of constraint */
      nvars = consdata->nvars;
      vars = consdata->vars;

      /* update number of SOS1 variables */
      for (i = 0; i < nvars; ++i)
      {
         SCIP_VAR* var;

         var = vars[i];

         if ( SCIPvarGetStatus(var) != SCIP_VARSTATUS_FIXED )
         {
            int ind;

            ind = SCIPvarGetIndex(var);
            assert( ind >= 0 && ind < ntotalvars );
            if ( ! nodecreated[ind] )
            {
               nodecreated[ind] = TRUE; /* mark node as counted */
               nodeorig[ind] = cntsos;
               ++cntsos;
            }
         }
      }
   }
   if ( cntsos <= 0 )
   {
      /* free buffer arrays */
      SCIPfreeBufferArray(scip, &nodecreated);
      SCIPfreeBufferArray(scip, &nodeorig);
      conshdlrdata->nsos1vars = 0;
      return SCIP_OKAY;
   }

   /* reinitialize vector 'nodecreated' */
   for (i = 0; i < ntotalvars; ++i)
      nodecreated[i] = FALSE;

   /* create conflict graph */
   SCIP_CALL( SCIPdigraphCreate(&conshdlrdata->conflictgraph, cntsos) );

   /* set up hash map */
   SCIP_CALL( SCIPhashmapCreate(&conshdlrdata->varhash, SCIPblkmem(scip), cntsos) );

   /* for every SOS1 constraint */
   cntsos = 0;
   for (c = 0; c < nconss; ++c)
   {
      SCIP_CONSDATA* consdata;
      SCIP_VAR** vars;
      int nvars;

      assert( conss[c] != NULL );

      /* get constraint data field of the constraint */
      consdata = SCIPconsGetData(conss[c]);
      assert( consdata != NULL );

      /* get variables and number of variables of constraint */
      nvars = consdata->nvars;
      vars = consdata->vars;

      /* add edges to the conflict graph and create node data for each of its nodes */
      for (i = 0; i < nvars; ++i)
      {
         SCIP_VAR* var;

         var = vars[i];

         if ( SCIPvarGetStatus(var) != SCIP_VARSTATUS_FIXED )
         {
            int indi;

            indi = SCIPvarGetIndex(var);

            if ( ! nodecreated[indi] )
            {
               SCIP_NODEDATA* nodedata = NULL;

               /* insert node number to hash map */
               assert( ! SCIPhashmapExists(conshdlrdata->varhash, var) );
               SCIP_CALL( SCIPhashmapInsert(conshdlrdata->varhash, var, (void*) (size_t) cntsos) );/*lint !e571*/
               assert( cntsos == (int) (size_t) SCIPhashmapGetImage(conshdlrdata->varhash, var) );
               assert( SCIPhashmapExists(conshdlrdata->varhash, var) );

               /* create node data */
               SCIP_CALL( SCIPallocMemory(scip, &nodedata) );
               nodedata->var = var;
               nodedata->lbboundvar = NULL;
               nodedata->ubboundvar = NULL;
               nodedata->lbboundcoef = 0.0;
               nodedata->ubboundcoef = 0.0;
               nodedata->lbboundcomp = FALSE;
               nodedata->ubboundcomp = FALSE;

               /* set node data */
               SCIPdigraphSetNodeData(conshdlrdata->conflictgraph, (void*)nodedata, cntsos);

               /* mark node and var data of node as created and update SOS1 counter */
               nodecreated[indi] = TRUE;
               ++cntsos;
            }

            /* add edges to the conflict graph */
            for (j = i+1; j < nvars; ++j)
            {
               var = vars[j];

               if ( SCIPvarGetStatus(var) != SCIP_VARSTATUS_FIXED )
               {
                  int indj;

                  indj = SCIPvarGetIndex(var);

                  /* in case indi = indj the variable will be deleted in the presolving step */
                  if ( indi != indj )
                  {
                     /* arcs have to be added 'safe' */
                     SCIP_CALL( SCIPdigraphAddArcSafe(conshdlrdata->conflictgraph, nodeorig[indi], nodeorig[indj], NULL) );
                     SCIP_CALL( SCIPdigraphAddArcSafe(conshdlrdata->conflictgraph, nodeorig[indj], nodeorig[indi], NULL) );
                  }
               }
            }
         }
      }
   }

   /* set number of problem variables that are contained in at least one SOS1 constraint */
   conshdlrdata->nsos1vars = cntsos;

   /* free buffer arrays */
   SCIPfreeBufferArray(scip, &nodecreated);
   SCIPfreeBufferArray(scip, &nodeorig);

   /* sort successors in ascending order */
   for (j = 0; j < conshdlrdata->nsos1vars; ++j)
   {
      int nsucc;

      nsucc = SCIPdigraphGetNSuccessors(conshdlrdata->conflictgraph, j);
      SCIPsortInt(SCIPdigraphGetSuccessors(conshdlrdata->conflictgraph, j), nsucc);
   }

   return SCIP_OKAY;
}


/** free conflict graph, nodedata and hashmap */
static
SCIP_RETCODE freeConflictgraph(
   SCIP_CONSHDLRDATA*    conshdlrdata        /**< constraint handler data */
   )
{
   int j;

   /* for every SOS1 variable */
   for (j = 0; j < conshdlrdata->nsos1vars; ++j)
   {
      SCIP_NODEDATA* nodedata;

      /* get node data */
      nodedata = (SCIP_NODEDATA*)SCIPdigraphGetNodeData(conshdlrdata->conflictgraph, j);
      assert( nodedata != NULL );

      /* free node data */
      SCIPfreeMemory(scip, &nodedata);
      SCIPdigraphSetNodeData(conshdlrdata->conflictgraph, NULL, j);
   }

   /* free conflict graph and hash map */
   if ( conshdlrdata->conflictgraph != NULL )
   {
      assert( conshdlrdata->nsos1vars > 0 );
      assert( conshdlrdata->varhash != NULL );
      SCIPhashmapFree(&conshdlrdata->varhash);
      SCIPdigraphFree(&conshdlrdata->conflictgraph);
   }

   return SCIP_OKAY;
}


/* ---------------------------- constraint handler callback methods ----------------------*/

/** copy method for constraint handler plugins (called when SCIP copies plugins) */
static
SCIP_DECL_CONSHDLRCOPY(conshdlrCopySOS1)
{  /*lint --e{715}*/
   assert( scip != NULL );
   assert( conshdlr != NULL );
   assert( strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0 );

   /* call inclusion method of constraint handler */
   SCIP_CALL( SCIPincludeConshdlrSOS1(scip) );

   *valid = TRUE;

   return SCIP_OKAY;
}


/** destructor of constraint handler to free constraint handler data (called when SCIP is exiting) */
static
SCIP_DECL_CONSFREE(consFreeSOS1)
{
   SCIP_CONSHDLRDATA* conshdlrdata;

   assert( scip != NULL );
   assert( conshdlr != NULL );
   assert( strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0 );

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   SCIPfreeMemory(scip, &conshdlrdata);

   return SCIP_OKAY;
}


/** solving process initialization method of constraint handler (called when branch and bound process is about to begin) */
static
SCIP_DECL_CONSINITSOL(consInitsolSOS1)
{  /*lint --e{715}*/
    SCIP_CONSHDLRDATA* conshdlrdata;

    assert( scip != NULL );
    assert( conshdlr != NULL );
    assert( strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0 );

    conshdlrdata = SCIPconshdlrGetData(conshdlr);
    assert( conshdlrdata != NULL );

    conshdlrdata->nsos1vars = 0;
    conshdlrdata->varhash = NULL;

    if ( nconss > 0 )
    {
       /* initialize conflict graph and hashmap for SOS1 variables */
       SCIP_CALL( initConflictgraph(scip, conshdlrdata, conss, nconss) );

       /* add data to conflict graph nodes */
       SCIP_CALL( setNodeDataSOS1(scip, conshdlr, conshdlrdata, nconss, conshdlrdata->nsos1vars) );
    }
    return SCIP_OKAY;
}


/** solving process deinitialization method of constraint handler (called before branch and bound process data is freed) */
static
SCIP_DECL_CONSEXITSOL(consExitsolSOS1)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;
   int c;

   assert( scip != NULL );
   assert( conshdlr != NULL );
   assert( strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0 );
   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert( conshdlrdata != NULL );

   /* free graph for storing local conflicts */
   if ( conshdlrdata->localconflicts != NULL )
      SCIPdigraphFree(&conshdlrdata->localconflicts);
   assert( conshdlrdata->localconflicts == NULL );

   /* free tclique graph and tclique data */
   if( conshdlrdata->tcliquegraph != NULL )
   {
      assert( conshdlrdata->tcliquedata != NULL );
      SCIPfreeMemory(scip, &conshdlrdata->tcliquedata);
      tcliqueFree(&conshdlrdata->tcliquegraph);
   }
   assert(conshdlrdata->tcliquegraph == NULL);
   assert(conshdlrdata->tcliquedata == NULL);

   /* free conflict graph */
   if ( nconss > 0 && conshdlrdata->nsos1vars > 0 )
   {
      SCIP_CALL( freeConflictgraph(conshdlrdata) );
   }
   assert( conshdlrdata->conflictgraph == NULL );

   /* check each constraint */
   for (c = 0; c < nconss; ++c)
   {
      SCIP_CONSDATA* consdata;

      assert( conss != NULL );
      assert( conss[c] != NULL );
      consdata = SCIPconsGetData(conss[c]);
      assert( consdata != NULL );

      SCIPdebugMessage("Exiting SOS1 constraint <%s>.\n", SCIPconsGetName(conss[c]) );

      /* free rows */
      if ( consdata->rowub != NULL )
         SCIP_CALL( SCIPreleaseRow(scip, &consdata->rowub) );

      if ( consdata->rowlb != NULL )
         SCIP_CALL( SCIPreleaseRow(scip, &consdata->rowlb) );
   }

   return SCIP_OKAY;
}


/** frees specific constraint data */
static
SCIP_DECL_CONSDELETE(consDeleteSOS1)
{
   assert( scip != NULL );
   assert( conshdlr != NULL );
   assert( cons != NULL );
   assert( consdata != NULL );
   assert( strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0 );

   SCIPdebugMessage("Deleting SOS1 constraint <%s>.\n", SCIPconsGetName(cons) );

   /* drop events on transformed variables */
   if ( SCIPconsIsTransformed(cons) )
   {
      SCIP_CONSHDLRDATA* conshdlrdata;
      int j;

      /* get constraint handler data */
      conshdlrdata = SCIPconshdlrGetData(conshdlr);
      assert( conshdlrdata != NULL );
      assert( conshdlrdata->eventhdlr != NULL );

      for (j = 0; j < (*consdata)->nvars; ++j)
      {
         SCIP_CALL( SCIPdropVarEvent(scip, (*consdata)->vars[j], SCIP_EVENTTYPE_BOUNDCHANGED, conshdlrdata->eventhdlr,
               (SCIP_EVENTDATA*)*consdata, -1) );
      }
   }

   SCIPfreeBlockMemoryArray(scip, &(*consdata)->vars, (*consdata)->maxvars);
   if ( (*consdata)->weights != NULL )
   {
      SCIPfreeBlockMemoryArray(scip, &(*consdata)->weights, (*consdata)->maxvars);
   }

   /* free rows */
   if ( (*consdata)->rowub != NULL )
   {
      SCIP_CALL( SCIPreleaseRow(scip, &(*consdata)->rowub) );
   }
   if ( (*consdata)->rowlb != NULL )
   {
      SCIP_CALL( SCIPreleaseRow(scip, &(*consdata)->rowlb) );
   }
   assert( (*consdata)->rowub == NULL );
   assert( (*consdata)->rowlb == NULL );

   SCIPfreeBlockMemory(scip, consdata);

   return SCIP_OKAY;
}


/** transforms constraint data into data belonging to the transformed problem */
static
SCIP_DECL_CONSTRANS(consTransSOS1)
{
   SCIP_CONSDATA* consdata;
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONSDATA* sourcedata;
   char s[SCIP_MAXSTRLEN];
   int j;

   assert( scip != NULL );
   assert( conshdlr != NULL );
   assert( strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0 );
   assert( sourcecons != NULL );
   assert( targetcons != NULL );

   /* get constraint handler data */
   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert( conshdlrdata != NULL );
   assert( conshdlrdata->eventhdlr != NULL );

   SCIPdebugMessage("Transforming SOS1 constraint: <%s>.\n", SCIPconsGetName(sourcecons) );

   /* get data of original constraint */
   sourcedata = SCIPconsGetData(sourcecons);
   assert( sourcedata != NULL );
   assert( sourcedata->nvars > 0 );
   assert( sourcedata->nvars <= sourcedata->maxvars );

   /* create constraint data */
   SCIP_CALL( SCIPallocBlockMemory(scip, &consdata) );

   consdata->nvars = sourcedata->nvars;
   consdata->maxvars = sourcedata->nvars;
   consdata->rowub = NULL;
   consdata->rowlb = NULL;
   consdata->nfixednonzeros = 0;
   consdata->local = sourcedata->local;
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &consdata->vars, consdata->nvars) );
   /* if weights were used */
   if ( sourcedata->weights != NULL )
   {
      SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &consdata->weights, sourcedata->weights, consdata->nvars) );
   }
   else
      consdata->weights = NULL;

   for (j = 0; j < sourcedata->nvars; ++j)
   {
      assert( sourcedata->vars[j] != 0 );
      SCIP_CALL( SCIPgetTransformedVar(scip, sourcedata->vars[j], &(consdata->vars[j])) );

      /* if variable is fixed to be nonzero */
      if ( SCIPisFeasPositive(scip, SCIPvarGetLbLocal(consdata->vars[j])) || SCIPisFeasNegative(scip, SCIPvarGetUbLocal(consdata->vars[j])) )
         ++(consdata->nfixednonzeros);
   }

   /* create transformed constraint with the same flags */
   (void) SCIPsnprintf(s, SCIP_MAXSTRLEN, "t_%s", SCIPconsGetName(sourcecons));
   SCIP_CALL( SCIPcreateCons(scip, targetcons, s, conshdlr, consdata,
         SCIPconsIsInitial(sourcecons), SCIPconsIsSeparated(sourcecons),
         SCIPconsIsEnforced(sourcecons), SCIPconsIsChecked(sourcecons),
         SCIPconsIsPropagated(sourcecons), SCIPconsIsLocal(sourcecons),
         SCIPconsIsModifiable(sourcecons), SCIPconsIsDynamic(sourcecons),
         SCIPconsIsRemovable(sourcecons), SCIPconsIsStickingAtNode(sourcecons)) );

   /* catch bound change events on variable */
   for (j = 0; j < consdata->nvars; ++j)
   {
      SCIP_CALL( SCIPcatchVarEvent(scip, consdata->vars[j], SCIP_EVENTTYPE_BOUNDCHANGED, conshdlrdata->eventhdlr,
            (SCIP_EVENTDATA*)consdata, NULL) );
   }

#ifdef SCIP_DEBUG
   if ( consdata->nfixednonzeros > 0 )
   {
      SCIPdebugMessage("constraint <%s> has %d variables fixed to be nonzero.\n", SCIPconsGetName(*targetcons),
         consdata->nfixednonzeros );
   }
#endif

   return SCIP_OKAY;
}


/** presolving method of constraint handler */
static
SCIP_DECL_CONSPRESOL(consPresolSOS1)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_EVENTHDLR* eventhdlr;
   SCIP_DIGRAPH* conflictgraph;
   SCIP_DIGRAPH* implgraph;
   SCIP_HASHMAP* implhash;
   SCIP_Bool* implnodes;
   SCIP_Bool** adjacencymatrix = NULL;
   SCIP_VAR** totalvars;
   int ntotalvars = 0;
   int nsos1vars;

   int oldnfixedvars;
   int oldndelconss;
   int oldnupgdconss;
   int nremovedvars;

   SCIP_Bool updateconfl;
   SCIP_Bool cutoff = FALSE;

   int* succ;
   int nsucc;
   int i;
   int j;

   assert( scip != NULL );
   assert( conshdlr != NULL );
   assert( strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0 );
   assert( result != NULL );

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert( conshdlrdata != NULL );

   SCIPdebugMessage("Presolving SOS1 constraints.\n");

   *result = SCIP_DIDNOTRUN;

   oldnfixedvars = *nfixedvars;
   oldndelconss = *ndelconss;
   oldnupgdconss = *nupgdconss;
   nremovedvars = 0;

   /* only run if success if possible */
   if( nconss > 0 && ( nrounds == 0 || nnewfixedvars > 0 || nnewaggrvars > 0 || nnewchgbds > 0 ) )
   {
      *result = SCIP_DIDNOTFIND;

      /* get constraint handler data */
      assert( SCIPconshdlrGetData(conshdlr) != NULL );
      eventhdlr = SCIPconshdlrGetData(conshdlr)->eventhdlr;
      assert( eventhdlr != NULL );

      /* initialize conflict graph */
      SCIP_CALL( initConflictgraph(scip, conshdlrdata, conss, nconss));

      /* get conflict graph */
      conflictgraph = conshdlrdata->conflictgraph;

      /* get number of SOS1 variables */
      nsos1vars = conshdlrdata->nsos1vars;

      if ( nsos1vars < 2 )
         return SCIP_OKAY;

      /* allocate buffer arrays */
      SCIP_CALL( SCIPallocBufferArray(scip, &adjacencymatrix, nsos1vars) );
      for (i = 0; i < nsos1vars; ++i)
         SCIP_CALL( SCIPallocBufferArray(scip, &adjacencymatrix[i], i+1) );

      /* create adjacency matrix */
      for (i = 0; i < nsos1vars; ++i)
      {
         for (j = 0; j < i+1; ++j)
            adjacencymatrix[i][j] = 0;
      }
      for (i = 0; i < nsos1vars; ++i)
      {
         succ = SCIPdigraphGetSuccessors(conflictgraph, i);
         nsucc = SCIPdigraphGetNSuccessors(conflictgraph, i);

         for (j = 0; j < nsucc; ++j)
         {
            if ( i > succ[j] )
               adjacencymatrix[i][succ[j]] = 1;
         }
      }

      /* perform one presolving round for SOS1 constraints */
      SCIP_CALL( presolRoundConssSOS1(scip, eventhdlr, conshdlr, conshdlrdata, conflictgraph, adjacencymatrix, conss, nconss, nsos1vars, naddconss, ndelconss, nupgdconss, nfixedvars, &nremovedvars, result) );

      /* if a cutoff has been detected */
      if ( *result == SCIP_CUTOFF )
      {
         /* free memory */
         for (j = nsos1vars-1; j >= 0; --j)
            SCIPfreeBufferArrayNull(scip, &adjacencymatrix[j]);
         SCIPfreeBufferArrayNull(scip, &adjacencymatrix);
         SCIP_CALL( freeConflictgraph(conshdlrdata));

         return SCIP_OKAY;
      }

      /* execute some pre-handling for bound tightening */
      if ( conshdlrdata->maxtightenbds != 0)
      {
         SCIP_VAR** probvars;
         int nprobvars;

         assert( ntotalvars == 0 );

         /* determine totalvars (union of SOS1 and problem variables) */
         probvars = SCIPgetVars(scip);
         nprobvars = SCIPgetNVars(scip);
         SCIP_CALL( SCIPhashmapCreate(&implhash, SCIPblkmem(scip), nsos1vars + nprobvars) );
         SCIP_CALL( SCIPallocBufferArray(scip, &totalvars, nsos1vars + nprobvars) );
         for (i = 0; i < nsos1vars; ++i)
         {
            SCIP_VAR* var;
            var = nodeGetVarSOS1(conflictgraph, i);

            /* insert node number to hash map */
            assert( ! SCIPhashmapExists(implhash, var) );
            SCIP_CALL( SCIPhashmapInsert(implhash, var, (void*) (size_t) ntotalvars) );/*lint !e571*/
            assert( ntotalvars == (int) (size_t) SCIPhashmapGetImage(implhash, var) );
            totalvars[ntotalvars++] = var;
         }
         for (i = 0; i < nprobvars; ++i)
         {
            SCIP_VAR* var;
            var = probvars[i];

            /* insert node number to hash map if not existent */
            if ( ! SCIPhashmapExists(implhash, var) )
            {
               SCIP_CALL( SCIPhashmapInsert(implhash, var, (void*) (size_t) ntotalvars) );/*lint !e571*/
               assert( ntotalvars == (int) (size_t) SCIPhashmapGetImage(implhash, var) );
               totalvars[ntotalvars++] = var;
            }
         }

         /* create implication graph */
         SCIP_CALL( SCIPdigraphCreate(&implgraph, ntotalvars) );
         SCIP_CALL( SCIPallocBufferArray(scip, &implnodes, nsos1vars) );
      }

      /* try to tighten the upper and lower bounds of the variables */
      updateconfl = FALSE;
      for (j = 0; (j < conshdlrdata->maxtightenbds || conshdlrdata->maxtightenbds == -1 ) && ! cutoff; ++j)
      {
         SCIP_Bool implupdate;
         int nchgbdssave;

         nchgbdssave = *nchgbds;

         assert( ntotalvars > 0 );
         SCIP_CALL( tightenVarsBoundsSOS1(scip, conshdlr, conflictgraph, implgraph, implhash, adjacencymatrix, totalvars, ntotalvars, nsos1vars, nchgbds, &implupdate, &cutoff) );
         if ( implupdate )
            updateconfl = TRUE;
         else if ( ! (*nchgbds > nchgbdssave ) )
            break;
      }

      if ( conshdlrdata->maxtightenbds != 0 )
      {
         /* try to update conflict graph based on the information of the implication graph */
         if ( updateconfl && conshdlrdata->updateconflpresol )
         {
            for (i = 0; i < nsos1vars; ++i)
            {
               for (j = 0; j < nsos1vars; ++j)
                  implnodes[j] = FALSE;
               SCIP_CALL( SCIPupdateConflictGraphSOS1(scip, conshdlr, conflictgraph, totalvars, implgraph, implhash, implnodes, adjacencymatrix, i, i, naddconss) );
            }
         }

         /* free memory needed for bound tightening */;
         SCIPhashmapFree(&implhash);
         SCIPfreeBufferArrayNull(scip, &totalvars);

         for (j = 0; j < ntotalvars; ++j)
         {
            SCIP_SUCCDATA** succdatas;
            int s;
            succdatas = (SCIP_SUCCDATA**) SCIPdigraphGetSuccessorsData(implgraph, j);
            nsucc = SCIPdigraphGetNSuccessors(implgraph, j);
            succ = SCIPdigraphGetSuccessors(implgraph, j);

            for (s = 0; s < nsucc; ++s)
               SCIPfreeMemory(scip, &succdatas[s]);
         }
         SCIPfreeBufferArrayNull(scip, &implnodes);
         SCIPdigraphFree(&implgraph);
      }

      /* free adjacency matrix */
      for (j = nsos1vars-1; j >= 0; --j)
         SCIPfreeBufferArrayNull(scip, &adjacencymatrix[j]);
      SCIPfreeBufferArrayNull(scip, &adjacencymatrix);

      /* if an infeasibility has been detected by bound tightening */
      if( cutoff )
      {
         SCIPdebugMessage("cutoff \n");
         *result = SCIP_CUTOFF;
      }

      /* free memory allocated in function initConflictgraph() */
      SCIP_CALL( freeConflictgraph(conshdlrdata));
   }
   (*nchgcoefs) += nremovedvars;

   SCIPdebugMessage("presolving fixed %d variables, removed %d variables, deleted %d constraints, and upgraded %d constraints.\n",
      *nfixedvars - oldnfixedvars, nremovedvars, *ndelconss - oldndelconss, *nupgdconss - oldnupgdconss);

   return SCIP_OKAY;
}


/** LP initialization method of constraint handler (called before the initial LP relaxation at a node is solved) */
static
SCIP_DECL_CONSINITLP(consInitlpSOS1)
{
   SCIP_CONSHDLRDATA* conshdlrdata;

   assert( scip != NULL );
   assert( conshdlr != NULL );
   assert( strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0 );

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert( conshdlrdata != NULL );

   /* checking for initial rows for SOS1 constraints */
   if( conshdlrdata->sepafromsos1 )
      SCIP_CALL( initsepaBoundInequalityFromSOS1Cons(scip, conshdlr, conshdlrdata, conss, nconss, NULL, FALSE, -1, NULL, NULL) );

   return SCIP_OKAY;
}


/** separation method of constraint handler for LP solutions */
static
SCIP_DECL_CONSSEPALP(consSepalpSOS1)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;
   int maxboundcuts;
   int ngen = 0;
   int depth;

   assert( scip != NULL );
   assert( conshdlr != NULL );
   assert( conss != NULL );
   assert( strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0 );
   assert( result != NULL );

   *result = SCIP_DIDNOTRUN;

   if ( nconss == 0 )
      return SCIP_OKAY;

   /* get constraint handler data */
   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert( conshdlrdata != NULL );

   /* check for boundcutsdepth < depth, maxboundcutsroot = 0 and maxboundcuts = 0 */
   depth = SCIPgetDepth(scip);
   if ( conshdlrdata->boundcutsdepth >= 0 && conshdlrdata->boundcutsdepth < depth )
      return SCIP_OKAY;

   /* only generate bound cuts if we are not close to terminating */
   if( SCIPisStopped(scip) )
      return SCIP_OKAY;

   /* determine maximal number of cuts*/
   if ( depth == 0 )
      maxboundcuts = conshdlrdata->maxboundcutsroot;
   else
      maxboundcuts = conshdlrdata->maxboundcuts;
   if ( maxboundcuts < 1 )
      return SCIP_OKAY;

   /* separate inequalities from SOS1 constraints */
   if( conshdlrdata->sepafromsos1 )
   {
      SCIP_CALL( initsepaBoundInequalityFromSOS1Cons(scip, conshdlr, conshdlrdata, conss, nconss, NULL, TRUE, maxboundcuts, &ngen, result) );
   }

   /* separate inequalities from the conflict graph */
   if( conshdlrdata->sepafromgraph )
   {
      SCIP_CALL( sepaBoundInequalitiesFromGraph(scip, conshdlr, conshdlrdata, NULL, maxboundcuts, &ngen, result) );
   }

   SCIPdebugMessage("Separated %d SOS1 constraints.\n", ngen);

   return SCIP_OKAY;
}


/** separation method of constraint handler for arbitrary primal solutions */
static
SCIP_DECL_CONSSEPASOL(consSepasolSOS1)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;
   int maxboundcuts;
   int ngen = 0;
   int depth;

   assert( scip != NULL );
   assert( conshdlr != NULL );
   assert( conss != NULL );
   assert( strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0 );
   assert( result != NULL );

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert( conshdlrdata != NULL );

   *result = SCIP_DIDNOTRUN;

   if ( nconss == 0 )
      return SCIP_OKAY;

   /* get constraint handler data */
   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert( conshdlrdata != NULL );

   /* check for boundcutsdepth < depth, maxboundcutsroot = 0 and maxboundcuts = 0 */
   depth = SCIPgetDepth(scip);
   if ( conshdlrdata->boundcutsdepth >= 0 && conshdlrdata->boundcutsdepth < depth )
      return SCIP_OKAY;

   /* only generate bound cuts if we are not close to terminating */
   if( SCIPisStopped(scip) )
      return SCIP_OKAY;

   /* determine maximal number of cuts*/
   if ( depth == 0 )
      maxboundcuts = conshdlrdata->maxboundcutsroot;
   else
      maxboundcuts = conshdlrdata->maxboundcuts;
   if ( maxboundcuts < 1 )
      return SCIP_OKAY;

   /* separate inequalities from sos1 constraints */
   if( conshdlrdata->sepafromsos1 )
      SCIP_CALL( initsepaBoundInequalityFromSOS1Cons(scip, conshdlr, conshdlrdata, conss, nconss, sol, TRUE, maxboundcuts, &ngen, result) );

   /* separate inequalities from the conflict graph */
   if( conshdlrdata->sepafromgraph )
   {
      SCIP_CALL( sepaBoundInequalitiesFromGraph(scip, conshdlr, conshdlrdata, sol, maxboundcuts, &ngen, result) );
   }

   SCIPdebugMessage("Separated %d SOS1 constraints.\n", ngen);

   return SCIP_OKAY;
}


/** constraint enforcing method of constraint handler for LP solutions */
static
SCIP_DECL_CONSENFOLP(consEnfolpSOS1)
{  /*lint --e{715}*/
   assert( scip != NULL );
   assert( conshdlr != NULL );
   assert( conss != NULL );
   assert( strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0 );
   assert( result != NULL );

   SCIP_CALL( enforceSOS1(scip, conshdlr, nconss, conss, result) );

   return SCIP_OKAY;
}


/** constraint enforcing method of constraint handler for pseudo solutions */
static
SCIP_DECL_CONSENFOPS(consEnfopsSOS1)
{  /*lint --e{715}*/
   assert( scip != NULL );
   assert( conshdlr != NULL );
   assert( conss != NULL );
   assert( strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0 );
   assert( result != NULL );

   SCIP_CALL( enforceSOS1(scip, conshdlr, nconss, conss, result) );

   return SCIP_OKAY;
}


/** feasibility check method of constraint handler for integral solutions
 *
 *  We simply check whether at most one variable is nonzero in the given solution.
 */
static
SCIP_DECL_CONSCHECK(consCheckSOS1)
{  /*lint --e{715}*/
   int c;

   assert( scip != NULL );
   assert( conshdlr != NULL );
   assert( conss != NULL );
   assert( strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0 );
   assert( result != NULL );

   /* check each constraint */
   for (c = 0; c < nconss; ++c)
   {
      SCIP_CONSDATA* consdata;
      int j;
      int cnt;

      cnt = 0;
      assert( conss[c] != NULL );
      consdata = SCIPconsGetData(conss[c]);
      assert( consdata != NULL );
      SCIPdebugMessage("Checking SOS1 constraint <%s>.\n", SCIPconsGetName(conss[c]));

      /* check all variables */
      for (j = 0; j < consdata->nvars; ++j)
      {
         /* if variable is nonzero */
         if ( ! SCIPisFeasZero(scip, SCIPgetSolVal(scip, sol, consdata->vars[j])) )
         {
            ++cnt;

            /* if more than one variable is nonzero */
            if ( cnt > 1 )
            {
               SCIP_CALL( SCIPresetConsAge(scip, conss[c]) );
               *result = SCIP_INFEASIBLE;

               if ( printreason )
               {
                  int l;

                  SCIP_CALL( SCIPprintCons(scip, conss[c], NULL) );
                  SCIPinfoMessage(scip, NULL, ";\nviolation: ");

                  for (l = 0; l < consdata->nvars; ++l)
                  {
                     /* if variable is nonzero */
                     if ( ! SCIPisFeasZero(scip, SCIPgetSolVal(scip, sol, consdata->vars[l])) )
                     {
                        SCIPinfoMessage(scip, NULL, "<%s> = %.15g ",
                           SCIPvarGetName(consdata->vars[l]), SCIPgetSolVal(scip, sol, consdata->vars[l]));
                     }
                  }
                  SCIPinfoMessage(scip, NULL, "\n");
               }
               return SCIP_OKAY;
            }
         }
      }
   }
   *result = SCIP_FEASIBLE;

   return SCIP_OKAY;
}


/** domain propagation method of constraint handler */
static
SCIP_DECL_CONSPROP(consPropSOS1)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_DIGRAPH* conflictgraph;
   int ngen = 0;

   assert( scip != NULL );
   assert( conshdlr != NULL );
   assert( conss != NULL );
   assert( strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0 );
   assert( result != NULL );
   *result = SCIP_DIDNOTRUN;

   assert( SCIPisTransformed(scip) );

   /* get constraint handler data */
   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert( conshdlrdata != NULL );

   /* get conflict graph */
   conflictgraph = conshdlrdata->conflictgraph;

   /* if conflict graph propagation shall be used */
   if ( conshdlrdata->conflictprop && conflictgraph != NULL )
   {
      int nsos1vars;
      int j;

      /* get number of SOS1 variables */
      nsos1vars = conshdlrdata->nsos1vars;

      /* check each SOS1 variable */
      for (j = 0; j < nsos1vars; ++j)
      {
         SCIP_VAR* var;

         var = nodeGetVarSOS1(conflictgraph, j);
         SCIPdebugMessage("Propagating SOS1 variable <%s>.\n", SCIPvarGetName(var) );

         /* if zero is outside the domain of variable */
         if ( SCIPisFeasPositive(scip, SCIPvarGetLbLocal(var)) || SCIPisFeasNegative(scip, SCIPvarGetUbLocal(var)) )
         {
            SCIP_VAR* succvar;
            SCIP_Real lb;
            SCIP_Real ub;
            int* succ;
            int nsucc;
            int s;

            /* fix all neighbors in the conflict graph to zero */
            succ = SCIPdigraphGetSuccessors(conflictgraph, j);
            nsucc = SCIPdigraphGetNSuccessors(conflictgraph, j);
            for (s = 0; s < nsucc; ++s)
            {
               succvar = nodeGetVarSOS1(conflictgraph, succ[s]);
               lb = SCIPvarGetLbLocal(succvar);
               ub = SCIPvarGetUbLocal(succvar);

               if ( ! SCIPisFeasZero(scip, lb) || ! SCIPisFeasZero(scip, ub) )
               {
                  /* if variable cannot be nonzero */
                  if ( SCIPisFeasPositive(scip, lb) || SCIPisFeasNegative(scip, ub) )
                  {
                     *result = SCIP_CUTOFF;
                     return SCIP_OKAY;
                  }

                  /* directly fix variable if it is not multi-aggregated */
                  if ( SCIPvarGetStatus(succvar) != SCIP_VARSTATUS_MULTAGGR )
                  {
                     SCIP_Bool infeasible;
                     SCIP_Bool tightened;

                     SCIP_CALL( SCIPtightenVarLb(scip, succvar, 0.0, FALSE, &infeasible, &tightened) );
                     assert( ! infeasible );
                     if ( tightened )
                        ++ngen;

                     SCIP_CALL( SCIPtightenVarUb(scip, succvar, 0.0, FALSE, &infeasible, &tightened) );
                     assert( ! infeasible );
                     if ( tightened )
                        ++ngen;
                  }
               }
            }
         }
      }
   }

   /* if SOS1 constraint propagation shall be used */
   if ( conshdlrdata->sosconsprop || conflictgraph == NULL )
   {
      int c;

      /* check each constraint */
      for (c = 0; c < nconss; ++c)
      {
         SCIP_CONS* cons;
         SCIP_CONSDATA* consdata;
         SCIP_Bool cutoff;

         *result = SCIP_DIDNOTFIND;
         assert( conss[c] != NULL );
         cons = conss[c];
         consdata = SCIPconsGetData(cons);
         assert( consdata != NULL );
         SCIPdebugMessage("Propagating SOS1 constraint <%s>.\n", SCIPconsGetName(cons) );

         *result = SCIP_DIDNOTFIND;
         SCIP_CALL( propSOS1(scip, cons, consdata, &cutoff, &ngen) );
         if ( cutoff )
         {
            *result = SCIP_CUTOFF;
            return SCIP_OKAY;
         }
      }
   }

   SCIPdebugMessage("Propagated %d domains.\n", ngen);
   if ( ngen > 0 )
      *result = SCIP_REDUCEDDOM;

   return SCIP_OKAY;
}


/** propagation conflict resolving method of constraint handler
 *
 *  We check which bound changes were the reason for infeasibility. We
 *  use that @a inferinfo stores the index of the variable that has
 *  bounds that fix it to be nonzero (these bounds are the reason). */
static
SCIP_DECL_CONSRESPROP(consRespropSOS1)
{  /*lint --e{715}*/
   SCIP_CONSDATA* consdata;
   SCIP_VAR* var;

   assert( scip != NULL );
   assert( cons != NULL );
   assert( strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0 );
   assert( infervar != NULL );
   assert( bdchgidx != NULL );
   assert( result != NULL );

   *result = SCIP_DIDNOTFIND;
   SCIPdebugMessage("Propagation resolution method of SOS1 constraint <%s>.\n", SCIPconsGetName(cons));

   consdata = SCIPconsGetData(cons);
   assert( consdata != NULL );
   assert( 0 <= inferinfo && inferinfo < consdata->nvars );
   var = consdata->vars[inferinfo];
   assert( var != infervar );

   /* check if lower bound of var was the reason */
   if ( SCIPisFeasPositive(scip, SCIPvarGetLbAtIndex(var, bdchgidx, FALSE)) )
   {
      SCIP_CALL( SCIPaddConflictLb(scip, var, bdchgidx) );
      *result = SCIP_SUCCESS;
   }

   /* check if upper bound of var was the reason */
   if ( SCIPisFeasNegative(scip, SCIPvarGetUbAtIndex(var, bdchgidx, FALSE)) )
   {
      SCIP_CALL( SCIPaddConflictUb(scip, var, bdchgidx) );
      *result = SCIP_SUCCESS;
   }

   return SCIP_OKAY;
}


/** variable rounding lock method of constraint handler
 *
 *  Let lb and ub be the lower and upper bounds of a
 *  variable. Preprocessing usually makes sure that lb <= 0 <= ub.
 *
 *  - If lb < 0 then rounding down may violate the constraint.
 *  - If ub > 0 then rounding up may violated the constraint.
 *  - If lb > 0 or ub < 0 then the constraint is infeasible and we do
 *    not have to deal with it here.
 *  - If lb == 0 then rounding down does not violate the constraint.
 *  - If ub == 0 then rounding up does not violate the constraint.
 */
static
SCIP_DECL_CONSLOCK(consLockSOS1)
{
   SCIP_CONSDATA* consdata;
   SCIP_VAR** vars;
   int nvars;
   int j;

   assert( scip != NULL );
   assert( conshdlr != NULL );
   assert( cons != NULL );
   assert( strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0 );
   consdata = SCIPconsGetData(cons);
   assert( consdata != NULL );

   SCIPdebugMessage("Locking constraint <%s>.\n", SCIPconsGetName(cons));

   vars = consdata->vars;
   nvars = consdata->nvars;
   assert( vars != NULL );

   for (j = 0; j < nvars; ++j)
   {
      SCIP_VAR* var;
      var = vars[j];

      /* if lower bound is negative, rounding down may violate constraint */
      if ( SCIPisFeasNegative(scip, SCIPvarGetLbLocal(var)) )
      {
         SCIP_CALL( SCIPaddVarLocks(scip, var, nlockspos, nlocksneg) );
      }

      /* additionally: if upper bound is positive, rounding up may violate constraint */
      if ( SCIPisFeasPositive(scip, SCIPvarGetUbLocal(var)) )
      {
         SCIP_CALL( SCIPaddVarLocks(scip, var, nlocksneg, nlockspos) );
      }
   }

   return SCIP_OKAY;
}


/** constraint display method of constraint handler */
static
SCIP_DECL_CONSPRINT(consPrintSOS1)
{  /*lint --e{715}*/
   SCIP_CONSDATA* consdata;
   int j;

   assert( scip != NULL );
   assert( conshdlr != NULL );
   assert( cons != NULL );
   assert( strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0 );

   consdata = SCIPconsGetData(cons);
   assert( consdata != NULL );

   for (j = 0; j < consdata->nvars; ++j)
   {
      if ( j > 0 )
         SCIPinfoMessage(scip, file, ", ");
      SCIP_CALL( SCIPwriteVarName(scip, file, consdata->vars[j], FALSE) );
      if ( consdata->weights == NULL )
         SCIPinfoMessage(scip, file, " (%d)", j+1);
      else
         SCIPinfoMessage(scip, file, " (%3.2f)", consdata->weights[j]);
   }

   return SCIP_OKAY;
}


/** constraint copying method of constraint handler */
static
SCIP_DECL_CONSCOPY(consCopySOS1)
{  /*lint --e{715}*/
   SCIP_CONSDATA* sourceconsdata;
   SCIP_VAR** sourcevars;
   SCIP_VAR** targetvars;
   SCIP_Real* sourceweights;
   SCIP_Real* targetweights;
   const char* consname;
   int nvars;
   int v;

   assert( scip != NULL );
   assert( sourcescip != NULL );
   assert( sourcecons != NULL );
   assert( strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(sourcecons)), CONSHDLR_NAME) == 0 );

   *valid = TRUE;

   if ( name != NULL )
      consname = name;
   else
      consname = SCIPconsGetName(sourcecons);

   SCIPdebugMessage("Copying SOS1 constraint <%s> ...\n", consname);

   sourceconsdata = SCIPconsGetData(sourcecons);
   assert( sourceconsdata != NULL );

   /* get variables and weights of the source constraint */
   nvars = sourceconsdata->nvars;

   if ( nvars == 0 )
      return SCIP_OKAY;

   sourcevars = sourceconsdata->vars;
   assert( sourcevars != NULL );
   sourceweights = sourceconsdata->weights;
   assert( sourceweights != NULL );

   /* duplicate variable array */
   SCIP_CALL( SCIPallocBufferArray(sourcescip, &targetvars, nvars) );
   SCIP_CALL( SCIPduplicateBufferArray(sourcescip, &targetweights, sourceweights, nvars) );

   /* get copied variables in target SCIP */
   for( v = 0; v < nvars && *valid; ++v )
   {
      SCIP_CALL( SCIPgetVarCopy(sourcescip, scip, sourcevars[v], &(targetvars[v]), varmap, consmap, global, valid) );
   }

    /* only create the target constraint, if all variables could be copied */
   if( *valid )
   {
      SCIP_CALL( SCIPcreateConsSOS1(scip, cons, consname, nvars, targetvars, targetweights,
            initial, separate, enforce, check, propagate, local, dynamic, removable, stickingatnode) );
   }

   /* free buffer array */
   SCIPfreeBufferArray(sourcescip, &targetweights);
   SCIPfreeBufferArray(sourcescip, &targetvars);

   return SCIP_OKAY;
}


/** constraint parsing method of constraint handler */
static
SCIP_DECL_CONSPARSE(consParseSOS1)
{  /*lint --e{715}*/
   SCIP_VAR* var;
   SCIP_Real weight;
   const char* s;
   char* t;

   *success = TRUE;
   s = str;

   /* create empty SOS1 constraint */
   SCIP_CALL( SCIPcreateConsSOS1(scip, cons, name, 0, NULL, NULL, initial, separate, enforce, check, propagate, local, dynamic, removable, stickingatnode) );

   /* loop through string */
   do
   {
      /* parse variable name */
      SCIP_CALL( SCIPparseVarName(scip, s, &var, &t) );
      s = t;

      /* skip until beginning of weight */
      while ( *s != '\0' && *s != '(' )
         ++s;

      if ( *s == '\0' )
      {
         SCIPverbMessage(scip, SCIP_VERBLEVEL_MINIMAL, NULL, "Syntax error: expected weight at input: %s\n", s);
         *success = FALSE;
         return SCIP_OKAY;
      }
      /* skip '(' */
      ++s;

      /* find weight */
      weight = strtod(s, &t);
      if ( t == NULL )
      {
         SCIPverbMessage(scip, SCIP_VERBLEVEL_MINIMAL, NULL, "Syntax error during parsing of the weight: %s\n", s);
         *success = FALSE;
         return SCIP_OKAY;
      }
      s = t;

      /* skip white space, ',', and ')' */
      while ( *s != '\0' && ( isspace((unsigned char)*s) ||  *s == ',' || *s == ')' ) )
         ++s;

      /* add variable */
      SCIP_CALL( SCIPaddVarSOS1(scip, *cons, var, weight) );
   }
   while ( *s != '\0' );

   return SCIP_OKAY;
}


/** constraint method of constraint handler which returns the variables (if possible) */
static
SCIP_DECL_CONSGETVARS(consGetVarsSOS1)
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
SCIP_DECL_CONSGETNVARS(consGetNVarsSOS1)
{  /*lint --e{715}*/
   SCIP_CONSDATA* consdata;

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   (*nvars) = consdata->nvars;
   (*success) = TRUE;

   return SCIP_OKAY;
}


/* ---------------- Callback methods of event handler ---------------- */

/* exec the event handler
 *
 * We update the number of variables fixed to be nonzero
 */
static
SCIP_DECL_EVENTEXEC(eventExecSOS1)
{
   SCIP_EVENTTYPE eventtype;
   SCIP_CONSDATA* consdata;
   SCIP_Real oldbound;
   SCIP_Real newbound;

   assert( eventhdlr != NULL );
   assert( eventdata != NULL );
   assert( strcmp(SCIPeventhdlrGetName(eventhdlr), EVENTHDLR_NAME) == 0 );
   assert( event != NULL );

   consdata = (SCIP_CONSDATA*)eventdata;
   assert( consdata != NULL );
   assert( 0 <= consdata->nfixednonzeros && consdata->nfixednonzeros <= consdata->nvars );

   oldbound = SCIPeventGetOldbound(event);
   newbound = SCIPeventGetNewbound(event);

   eventtype = SCIPeventGetType(event);
   switch ( eventtype )
   {
   case SCIP_EVENTTYPE_LBTIGHTENED:
      /* if variable is now fixed to be nonzero */
      if ( ! SCIPisFeasPositive(scip, oldbound) && SCIPisFeasPositive(scip, newbound) )
         ++(consdata->nfixednonzeros);
      break;
   case SCIP_EVENTTYPE_UBTIGHTENED:
      /* if variable is now fixed to be nonzero */
      if ( ! SCIPisFeasNegative(scip, oldbound) && SCIPisFeasNegative(scip, newbound) )
         ++(consdata->nfixednonzeros);
      break;
   case SCIP_EVENTTYPE_LBRELAXED:
      /* if variable is not fixed to be nonzero anymore */
      if ( SCIPisFeasPositive(scip, oldbound) && ! SCIPisFeasPositive(scip, newbound) )
         --(consdata->nfixednonzeros);
      break;
   case SCIP_EVENTTYPE_UBRELAXED:
      /* if variable is not fixed to be nonzero anymore */
      if ( SCIPisFeasNegative(scip, oldbound) && ! SCIPisFeasNegative(scip, newbound) )
         --(consdata->nfixednonzeros);
      break;
   default:
      SCIPerrorMessage("invalid event type.\n");
      return SCIP_INVALIDDATA;
   }
   assert( 0 <= consdata->nfixednonzeros && consdata->nfixednonzeros <= consdata->nvars );

   SCIPdebugMessage("changed bound of variable <%s> from %f to %f (nfixednonzeros: %d).\n", SCIPvarGetName(SCIPeventGetVar(event)),
                    oldbound, newbound, consdata->nfixednonzeros);

   return SCIP_OKAY;
}


/* ---------------- Constraint specific interface methods ---------------- */

/** creates the handler for SOS1 constraints and includes it in SCIP */
SCIP_RETCODE SCIPincludeConshdlrSOS1(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONSHDLR* conshdlr;

   /* create constraint handler data */
   SCIP_CALL( SCIPallocMemory(scip, &conshdlrdata) );
   conshdlrdata->branchsos = TRUE;
   conshdlrdata->eventhdlr = NULL;
   conshdlrdata->conflictgraph = NULL;
   conshdlrdata->localconflicts = NULL;
   conshdlrdata->isconflocal = FALSE;
   conshdlrdata->nboundcuts = 0;
   conshdlrdata->tcliquegraph = NULL;
   conshdlrdata->tcliquedata = NULL;

   /* create event handler for bound change events */
   SCIP_CALL( SCIPincludeEventhdlrBasic(scip, &conshdlrdata->eventhdlr, EVENTHDLR_NAME, EVENTHDLR_DESC, eventExecSOS1, NULL) );
   if ( conshdlrdata->eventhdlr == NULL )
   {
      SCIPerrorMessage("event handler for SOS1 constraints not found.\n");
      return SCIP_PLUGINNOTFOUND;
   }

   /* include constraint handler */
   SCIP_CALL( SCIPincludeConshdlrBasic(scip, &conshdlr, CONSHDLR_NAME, CONSHDLR_DESC,
         CONSHDLR_ENFOPRIORITY, CONSHDLR_CHECKPRIORITY, CONSHDLR_EAGERFREQ, CONSHDLR_NEEDSCONS,
         consEnfolpSOS1, consEnfopsSOS1, consCheckSOS1, consLockSOS1, conshdlrdata) );
   assert(conshdlr != NULL);

   /* set non-fundamental callbacks via specific setter functions */
   SCIP_CALL( SCIPsetConshdlrCopy(scip, conshdlr, conshdlrCopySOS1, consCopySOS1) );
   SCIP_CALL( SCIPsetConshdlrDelete(scip, conshdlr, consDeleteSOS1) );
   SCIP_CALL( SCIPsetConshdlrExitsol(scip, conshdlr, consExitsolSOS1) );
   SCIP_CALL( SCIPsetConshdlrInitsol(scip, conshdlr, consInitsolSOS1) );
   SCIP_CALL( SCIPsetConshdlrFree(scip, conshdlr, consFreeSOS1) );
   SCIP_CALL( SCIPsetConshdlrGetVars(scip, conshdlr, consGetVarsSOS1) );
   SCIP_CALL( SCIPsetConshdlrGetNVars(scip, conshdlr, consGetNVarsSOS1) );
   SCIP_CALL( SCIPsetConshdlrInitlp(scip, conshdlr, consInitlpSOS1) );
   SCIP_CALL( SCIPsetConshdlrParse(scip, conshdlr, consParseSOS1) );
   SCIP_CALL( SCIPsetConshdlrPresol(scip, conshdlr, consPresolSOS1, CONSHDLR_MAXPREROUNDS, CONSHDLR_DELAYPRESOL) );
   SCIP_CALL( SCIPsetConshdlrPrint(scip, conshdlr, consPrintSOS1) );
   SCIP_CALL( SCIPsetConshdlrProp(scip, conshdlr, consPropSOS1, CONSHDLR_PROPFREQ, CONSHDLR_DELAYPROP, CONSHDLR_PROP_TIMING) );
   SCIP_CALL( SCIPsetConshdlrResprop(scip, conshdlr, consRespropSOS1) );
   SCIP_CALL( SCIPsetConshdlrSepa(scip, conshdlr, consSepalpSOS1, consSepasolSOS1, CONSHDLR_SEPAFREQ, CONSHDLR_SEPAPRIORITY, CONSHDLR_DELAYSEPA) );
   SCIP_CALL( SCIPsetConshdlrTrans(scip, conshdlr, consTransSOS1) );

   /* add SOS1 constraint handler parameters */

   /* presolving parameters */
   SCIP_CALL( SCIPaddIntParam(scip, "constraints/"CONSHDLR_NAME"/maxconsdelayext",
         "delay clique extension if number of sos1 constraints is larger than predefined value (-1: no limit)",
         &conshdlrdata->maxconsdelayext, TRUE, DEFAULT_MAXCONSDELAYEXT, -1, INT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(scip, "constraints/"CONSHDLR_NAME"/maxextensions",
         "maximal number of extensions that will be computed for each SOS1 constraint  (-1: no limit)",
         &conshdlrdata->maxextensions, TRUE, DEFAULT_MAXEXTENSIONS, -1, INT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(scip, "constraints/"CONSHDLR_NAME"/maxtightenbds",
         "maximal number of bound tightening rounds per presolving round (-1: no limit)",
         &conshdlrdata->maxtightenbds, TRUE, DEFAULT_MAXTIGHTENBDS, -1, INT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "constraints/"CONSHDLR_NAME"/updateconflpresol",
         "if TRUE then update conflict graph during presolving procedure",
         &conshdlrdata->updateconflpresol, TRUE, DEFAULT_UPDATECONFLPRESOL, NULL, NULL) );

   /* propagation parameters */
   SCIP_CALL( SCIPaddBoolParam(scip, "constraints/"CONSHDLR_NAME"/conflictprop",
         "whether to use conflict graph propagation",
         &conshdlrdata->conflictprop, TRUE, DEFAULT_CONFLICTPROP, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "constraints/"CONSHDLR_NAME"/sosconsprop",
         "whether to use SOS1 constraint propagation",
         &conshdlrdata->sosconsprop, TRUE, DEFAULT_SOSCONSPROP, NULL, NULL) );

   /* branching parameters */
   SCIP_CALL( SCIPaddBoolParam(scip, "constraints/"CONSHDLR_NAME"/branchsos",
         "Use SOS1 branching in enforcing (otherwise leave decision to branching rules)?",
         &conshdlrdata->branchsos, FALSE, TRUE, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "constraints/"CONSHDLR_NAME"/branchnonzeros",
         "Branch on SOS constraint with most number of nonzeros?",
         &conshdlrdata->branchnonzeros, FALSE, FALSE, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "constraints/"CONSHDLR_NAME"/branchweight",
         "Branch on SOS cons. with highest nonzero-variable weight for branching (needs branchnonzeros = false)?",
         &conshdlrdata->branchweight, FALSE, FALSE, NULL, NULL) );

   /* separation parameters */
   SCIP_CALL( SCIPaddBoolParam(scip, "constraints/"CONSHDLR_NAME"/sepafromsos1",
         "if TRUE separate bound inequalities from initial SOS1 constraints",
         &conshdlrdata->sepafromsos1, TRUE, DEFAULT_SEPAFROMSOS1, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "constraints/"CONSHDLR_NAME"/sepafromgraph",
         "if TRUE separate bound inequalities from the conflict graph",
         &conshdlrdata->sepafromgraph, TRUE, DEFAULT_SEPAFROMGRAPH, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(scip, "constraints/"CONSHDLR_NAME"/boundcutsdepth",
         "node depth of separating bound cuts (-1: no limit)",
         &conshdlrdata->boundcutsdepth, TRUE, DEFAULT_BOUNDCUTSDEPTH, -1, INT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(scip, "constraints/"CONSHDLR_NAME"/maxboundcuts",
         "maximal number of bound cuts separated per branching node",
         &conshdlrdata->maxboundcuts, TRUE, DEFAULT_MAXBOUNDCUTS, 0, INT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(scip, "constraints/"CONSHDLR_NAME"/maxboundcutsroot",
         "maximal number of bound cuts separated per iteration in the root node",
         &conshdlrdata->maxboundcutsroot, TRUE, DEFAULT_MAXBOUNDCUTSROOT, 0, INT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "constraints/"CONSHDLR_NAME"/strthenboundcuts",
         "if TRUE then bound cuts are strengthened in case bound variables are available",
         &conshdlrdata->strthenboundcuts, TRUE, DEFAULT_STRTHENBOUNDCUTS, NULL, NULL) );

   return SCIP_OKAY;
}


/** creates and captures a SOS1 constraint
 *
 *  We set the constraint to not be modifable. If the weights are non
 *  NULL, the variables are ordered according to these weights (in
 *  ascending order).
 *
 *  @note the constraint gets captured, hence at one point you have to release it using the method SCIPreleaseCons()
 */
SCIP_RETCODE SCIPcreateConsSOS1(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS**           cons,               /**< pointer to hold the created constraint */
   const char*           name,               /**< name of constraint */
   int                   nvars,              /**< number of variables in the constraint */
   SCIP_VAR**            vars,               /**< array with variables of constraint entries */
   SCIP_Real*            weights,            /**< weights determining the variable order, or NULL if natural order should be used */
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
   SCIP_Bool modifiable;
   SCIP_Bool transformed;
   int v;

   modifiable = FALSE;

   /* find the SOS1 constraint handler */
   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);
   if ( conshdlr == NULL )
   {
      SCIPerrorMessage("<%s> constraint handler not found\n", CONSHDLR_NAME);
      return SCIP_PLUGINNOTFOUND;
   }

   /* are we in the transformed problem? */
   transformed = SCIPgetStage(scip) >= SCIP_STAGE_TRANSFORMED;

   /* create constraint data */
   SCIP_CALL( SCIPallocBlockMemory(scip, &consdata) );
   consdata->vars = NULL;
   consdata->nvars = nvars;
   consdata->maxvars = nvars;
   consdata->rowub = NULL;
   consdata->rowlb = NULL;
   consdata->nfixednonzeros = transformed ? 0 : -1;
   consdata->weights = NULL;
   consdata->local = local;

   if ( nvars > 0 )
   {
      SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &consdata->vars, vars, nvars) );

      /* check weights */
      if ( weights != NULL )
      {
         /* store weights */
         SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &consdata->weights, weights, nvars) );

         /* sort variables - ascending order */
         SCIPsortRealPtr(consdata->weights, (void**)consdata->vars, nvars);
      }
   }
   else
   {
      assert( weights == NULL );
   }

   /* branching on multiaggregated variables does not seem to work well, so avoid it */
   for (v = 0; v < nvars; ++v)
      SCIP_CALL( SCIPmarkDoNotMultaggrVar(scip, consdata->vars[v]) );

   /* create constraint */
   SCIP_CALL( SCIPcreateCons(scip, cons, name, conshdlr, consdata, initial, separate, enforce, check, propagate,
         local, modifiable, dynamic, removable, stickingatnode) );
   assert(transformed == SCIPconsIsTransformed(*cons));

   /* replace original variables by transformed variables in transformed constraint, add locks, and catch events */
   for( v = nvars - 1; v >= 0; --v )
   {
      /* always use transformed variables in transformed constraints */
      if ( transformed )
      {
         SCIP_CALL( SCIPgetTransformedVar(scip, consdata->vars[v], &(consdata->vars[v])) );
      }
      assert( consdata->vars[v] != NULL );
      assert( transformed == SCIPvarIsTransformed(consdata->vars[v]) );

      /* handle the new variable */
      SCIP_CALL( handleNewVariableSOS1(scip, *cons, consdata, consdata->vars[v], transformed) );
   }

   return SCIP_OKAY;
}


/** creates and captures a SOS1 constraint with all constraint flags set to their default values.
 *
 *  @warning Do NOT set the constraint to be modifiable manually, because this might lead
 *  to wrong results as the variable array will not be resorted
 *
 *  @note the constraint gets captured, hence at one point you have to release it using the method SCIPreleaseCons()
 */
SCIP_RETCODE SCIPcreateConsBasicSOS1(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS**           cons,               /**< pointer to hold the created constraint */
   const char*           name,               /**< name of constraint */
   int                   nvars,              /**< number of variables in the constraint */
   SCIP_VAR**            vars,               /**< array with variables of constraint entries */
   SCIP_Real*            weights             /**< weights determining the variable order, or NULL if natural order should be used */
   )
{
   SCIP_CALL( SCIPcreateConsSOS1( scip, cons, name, nvars, vars, weights, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   return SCIP_OKAY;
}


/** adds variable to SOS1 constraint, the position is determined by the given weight */
SCIP_RETCODE SCIPaddVarSOS1(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint */
   SCIP_VAR*             var,                /**< variable to add to the constraint */
   SCIP_Real             weight              /**< weight determining position of variable */
   )
{
   assert( scip != NULL );
   assert( var != NULL );
   assert( cons != NULL );

   SCIPdebugMessage("adding variable <%s> to constraint <%s> with weight %g\n", SCIPvarGetName(var), SCIPconsGetName(cons), weight);

   if ( strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(cons)), CONSHDLR_NAME) != 0 )
   {
      SCIPerrorMessage("constraint is not an SOS1 constraint.\n");
      return SCIP_INVALIDDATA;
   }

   SCIP_CALL( addVarSOS1(scip, cons, var, weight) );

   return SCIP_OKAY;
}


/** appends variable to SOS1 constraint */
SCIP_RETCODE SCIPappendVarSOS1(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint */
   SCIP_VAR*             var                 /**< variable to add to the constraint */
   )
{
   assert( scip != NULL );
   assert( var != NULL );
   assert( cons != NULL );

   SCIPdebugMessage("appending variable <%s> to constraint <%s>\n", SCIPvarGetName(var), SCIPconsGetName(cons));

   if ( strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(cons)), CONSHDLR_NAME) != 0 )
   {
      SCIPerrorMessage("constraint is not an SOS1 constraint.\n");
      return SCIP_INVALIDDATA;
   }

   SCIP_CALL( appendVarSOS1(scip, cons, var) );

   return SCIP_OKAY;
}


/** gets number of variables in SOS1 constraint */
int SCIPgetNVarsSOS1(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint */
   )
{
   SCIP_CONSDATA* consdata;

   assert( scip != NULL );
   assert( cons != NULL );

   if ( strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(cons)), CONSHDLR_NAME) != 0 )
   {
      SCIPerrorMessage("constraint is not an SOS1 constraint.\n");
      SCIPABORT();
      return -1;  /*lint !e527*/
   }

   consdata = SCIPconsGetData(cons);
   assert( consdata != NULL );

   return consdata->nvars;
}


/** gets array of variables in SOS1 constraint */
SCIP_VAR** SCIPgetVarsSOS1(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint data */
   )
{
   SCIP_CONSDATA* consdata;

   assert( scip != NULL );
   assert( cons != NULL );

   if ( strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(cons)), CONSHDLR_NAME) != 0 )
   {
      SCIPerrorMessage("constraint is not an SOS1 constraint.\n");
      SCIPABORT();
      return NULL;  /*lint !e527*/
   }

   consdata = SCIPconsGetData(cons);
   assert( consdata != NULL );

   return consdata->vars;
}


/** gets array of weights in SOS1 constraint (or NULL if not existent) */
SCIP_Real* SCIPgetWeightsSOS1(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint data */
   )
{
   SCIP_CONSDATA* consdata;

   assert( scip != NULL );
   assert( cons != NULL );

   if ( strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(cons)), CONSHDLR_NAME) != 0 )
   {
      SCIPerrorMessage("constraint is not an SOS1 constraint.\n");
      SCIPABORT();
      return NULL;  /*lint !e527*/
   }

   consdata = SCIPconsGetData(cons);
   assert( consdata != NULL );

   return consdata->weights;
}


/** gets conflict graph of SOS1 constraints (or NULL if not existent)
 *
 *  Note: The conflict graph is globally valid; local changes are not taken into account.
 */
SCIP_DIGRAPH* SCIPgetConflictgraphSOS1(
   SCIP_CONSHDLR*        conshdlr            /**< SOS1 constraint handler */
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata;

   assert( conshdlr != NULL );

   if ( strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) != 0 )
   {
      SCIPerrorMessage("not an SOS1 constraint handler.\n");
      SCIPABORT();
   }
   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert( conshdlrdata != NULL );

   return conshdlrdata->conflictgraph;
}


/** gets number of problem variables that are involved in at least one SOS1 constraint */
int SCIPgetNSOS1Vars(
   SCIP_CONSHDLR*        conshdlr            /**< SOS1 constraint handler */
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata;

   assert( conshdlr != NULL );

   if ( strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) != 0 )
   {
      SCIPerrorMessage("not an SOS1 constraint handler.\n");
      SCIPABORT();
   }
   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert( conshdlrdata != NULL );

   return conshdlrdata->nsos1vars;
}


/** returns whether variable is involved in an SOS1 constraint */
SCIP_Bool varIsSOS1(
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler */
   SCIP_VAR*             var                 /**< variable */
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata;

   assert( var != NULL );
   assert( conshdlr != NULL );

   if ( strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) != 0 )
   {
      SCIPerrorMessage("not an SOS1 constraint handler.\n");
      SCIPABORT();
   }
   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert( conshdlrdata != NULL );

   if ( conshdlrdata->varhash == NULL || ! SCIPhashmapExists(conshdlrdata->varhash, var) )
      return FALSE;

   return TRUE;
}


/** returns SOS1 index of variable or -1 if variable is not involved in an SOS1 constraint */
int varGetNodeSOS1(
   SCIP_CONSHDLR*        conshdlr,            /**< SOS1 constraint handler */
   SCIP_VAR*             var                  /**< variable */
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata;

   assert( conshdlr != NULL );
   assert( var != NULL );

   if ( strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) != 0 )
   {
      SCIPerrorMessage("not an SOS1 constraint handler.\n");
      SCIPABORT();
   }
   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert( conshdlrdata != NULL );

   if ( ! SCIPhashmapExists(conshdlrdata->varhash, var) )
      return -1;

   return (int) (size_t) SCIPhashmapGetImage(conshdlrdata->varhash, var);
}


/** returns variable that belongs to a given node from the conflictgraph */
SCIP_VAR* nodeGetVarSOS1(
   SCIP_DIGRAPH*         conflictgraph,      /**< conflict graph */
   int                   node                /**< node from the conflict graph */
   )
{
   SCIP_NODEDATA* nodedata;

   assert( conflictgraph != NULL );
   assert( node >= 0 && node < SCIPdigraphGetNNodes(conflictgraph) );

   /* get node data */
   nodedata = (SCIP_NODEDATA*)SCIPdigraphGetNodeData(conflictgraph, node);

   if ( nodedata == NULL )
   {
      SCIPerrorMessage("variable is not assigned to an index.\n");
      SCIPABORT();
   }

   return nodedata->var;
}


/** returns solution value of imaginary binary big-M variable of a given node from the conflict graph */
SCIP_Real SCIPnodeGetSolvalBinaryBigMSOS1(
   SCIP*                 scip,               /**< SCIP pointer */
   SCIP_DIGRAPH*         conflictgraph,      /**< conflict graph */
   SCIP_SOL*         	 sol,                /**< primal solution, or NULL for current LP/pseudo solution */
   int                   node                /**< node of the conflict graph */
   )
{
   SCIP_Real bound;
   SCIP_VAR* var;
   SCIP_Real val;

   assert( scip != NULL );
   assert( conflictgraph != NULL );
   assert( node >= 0 && node < SCIPdigraphGetNNodes(conflictgraph) );

   var = nodeGetVarSOS1(conflictgraph, node);
   val = SCIPgetSolVal(scip, sol, var);

   if ( SCIPisFeasNegative(scip, val) )
   {
      bound = SCIPvarGetLbLocal(var);
      assert( SCIPisFeasNegative(scip, bound) );

      if ( SCIPisInfinity(scip, -val) )
         return 1.0;
      else if ( SCIPisInfinity(scip, -bound) )
         return 0.0;
      else
         return (val/bound);
   }
   else if ( SCIPisFeasPositive(scip, val) )
   {
      bound = SCIPvarGetUbLocal(var);
      assert( SCIPisFeasPositive(scip, bound) );
      assert( SCIPisFeasPositive(scip, val) );

      if ( SCIPisInfinity(scip, val) )
         return 1.0;
      else if ( SCIPisInfinity(scip, bound) )
         return 0.0;
      else
         return (val/bound);
   }
   else
      return 0.0;
}


/** gets (variable) lower bound value of current LP relaxation solution for a given node from the conflict graph */
SCIP_Real SCIPnodeGetSolvalVarboundLbSOS1(
   SCIP*                 scip,               /**< SCIP pointer */
   SCIP_DIGRAPH*         conflictgraph,      /**< conflict graph */
   SCIP_SOL*         	 sol,                /**< primal solution, or NULL for current LP/pseudo solution */
   int                   node                /**< node of the conflict graph */
   )
{
   SCIP_NODEDATA* nodedata;

   assert( scip != NULL );
   assert( conflictgraph != NULL );
   assert( node >= 0 && node < SCIPdigraphGetNNodes(conflictgraph) );

   /* get node data */
   nodedata = (SCIP_NODEDATA*)SCIPdigraphGetNodeData(conflictgraph, node);
   assert( nodedata != NULL );

   /* if variable is not involved in a variable upper bound constraint */
   if ( nodedata->lbboundvar == NULL || ! nodedata->lbboundcomp )
      return SCIPvarGetLbLocal(nodedata->var);

   return nodedata->lbboundcoef * SCIPgetSolVal(scip, sol, nodedata->lbboundvar);
}


/** gets (variable) upper bound value of current LP relaxation solution for a given node from the conflict graph */
SCIP_Real SCIPnodeGetSolvalVarboundUbSOS1(
   SCIP*                 scip,               /**< SCIP pointer */
   SCIP_DIGRAPH*         conflictgraph,      /**< conflict graph */
   SCIP_SOL*         	 sol,                /**< primal solution, or NULL for current LP/pseudo solution */
   int                   node                /**< node of the conflict graph */
   )
{
   SCIP_NODEDATA* nodedata;

   assert( scip != NULL );
   assert( conflictgraph != NULL );
   assert( node >= 0 && node < SCIPdigraphGetNNodes(conflictgraph) );

   /* get node data */
   nodedata = (SCIP_NODEDATA*)SCIPdigraphGetNodeData(conflictgraph, node);
   assert( nodedata != NULL );

   /* if variable is not involved in a variable upper bound constraint */
   if ( nodedata->ubboundvar == NULL || ! nodedata->ubboundcomp )
      return SCIPvarGetUbLocal(nodedata->var);

   return nodedata->ubboundcoef * SCIPgetSolVal(scip, sol, nodedata->ubboundvar);
}
