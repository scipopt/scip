/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2019 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not visit scip.zib.de.         */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   presol_milp.cpp
 * @brief  MILP presolver
 * @author Leona Gottwald
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "scip/presol_milp.h"

#ifndef SCIP_WITH_PRESOLVELIB

/** creates the xyz presolver and includes it in SCIP */
SCIP_RETCODE SCIPincludePresolMILP(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   return SCIP_OKAY;
}

#else

#include <assert.h>
#include "scip/cons_linear.h"
#include "scip/pub_matrix.h"
#include "scip/pub_presol.h"
#include "scip/pub_var.h"
#include "scip/pub_cons.h"
#include "scip/pub_message.h"
#include "scip/scip_presol.h"
#include "scip/scip_var.h"
#include "scip/scip_mem.h"
#include "scip/scip_prob.h"
#include "scip/scip_param.h"
#include "scip/scip_cons.h"
#include "scip/scip_numerics.h"
#include "scip/scip_timing.h"
#include "scip/scip_message.h"
#include "core/Presolve.hpp"
#include "core/ProblemBuilder.hpp"
#include "tbb/task_scheduler_init.h"


#define PRESOL_NAME            "milp"
#define PRESOL_DESC            "MILP specific presolving methods"
#define PRESOL_PRIORITY        -9999999 /**< priority of the presolver (>= 0: before, < 0: after constraint handlers); combined with propagators */
#define PRESOL_MAXROUNDS             -1 /**< maximal number of presolving rounds the presolver participates in (-1: no limit) */
#define PRESOL_TIMING           SCIP_PRESOLTIMING_MEDIUM /* timing of the presolver (fast, medium, or exhaustive) */

/* default parameter values */
#define DEFAULT_THREADS            1         /**< maximum number of threads presolving may use (0: automatic) */
#define DEFAULT_MAXFILLINPERSUBST  8         /**< maximal possible fillin for substitutions to be considered */
#define DEFAULT_MODIFYCONSFAC      0.0       /**< modify SCIP constraints when the number of nonzeros is at most this
                                              *   factor times the number of nonzeros before presolving */
#define DEFAULT_ENABLEPARALLELROWS TRUE      /**< should the parallel rows presolver be enabled within the presolve library? */
#define DEFAULT_ENABLEDOMCOL       TRUE      /**< should the dominated column presolver be enabled within the presolve library? */
#define DEFAULT_ENABLEDUALINFER    TRUE      /**< should the dualinfer presolver be enabled within the presolve library? */
#define DEFAULT_ENABLEMULTIAGGR    TRUE      /**< should the multi-aggregation presolver be enabled within the presolve library? */
#define DEFAULT_ENABLEPROBING      TRUE      /**< should the probing presolver be enabled within the presolve library? */
#define DEFAULT_ENABLESPARSIFY     FALSE     /**< should the sparsify presolver be enabled within the presolve library? */

/*
 * Data structures
 */

/** presolver data */
struct SCIP_PresolData
{
   tbb::task_scheduler_init schedulerinit;   /**< initialization object for tbb scheduler */
   int lastncols;                            /**< the number of columns from the last call */
   int lastnrows;                            /**< the number of rows from the last call */
   int threads;                              /**< maximum number of threads presolving may use (0: automatic) */
   int maxfillinpersubstitution;             /**< maximal possible fillin for substitutions to be considered */
   SCIP_Bool enablesparsify;                 /**< should the sparsify presolver be enabled within the presolve library? */
   SCIP_Bool enabledomcol;                   /**< should the dominated column presolver be enabled within the presolve library? */
   SCIP_Bool enableprobing;                  /**< should the probing presolver be enabled within the presolve library? */
   SCIP_Bool enabledualinfer;                /**< should the dualinfer presolver be enabled within the presolve library? */
   SCIP_Bool enablemultiaggr;                /**< should the multi-aggregation presolver be enabled within the presolve library? */
   SCIP_Bool enableparallelrows;             /**< should the parallel rows presolver be enabled within the presolve library? */
   SCIP_Real modifyconsfac;                  /**< modify SCIP constraints when the number of nonzeros is at most this
                                              *   factor times the number of nonzeros before presolving */
};


/*
 * Local methods
 */

/** builds the presolvelib problem datastructure from the matrix */
static
Problem<SCIP_Real> buildProblem(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_MATRIX*          matrix              /**< initialized SCIP_MATRIX data structure */
   )
{
   ProblemBuilder<SCIP_Real> builder;

   /* build problem from matrix */
   int nnz = SCIPmatrixGetNNonzs(matrix);
   int ncols = SCIPmatrixGetNColumns(matrix);
   int nrows = SCIPmatrixGetNRows(matrix);
   builder.reserve(nnz, nrows, ncols);

   /* set up columns */
   builder.setNumCols(ncols);
   for(int i = 0; i != ncols; ++i)
   {
      SCIP_VAR* var = SCIPmatrixGetVar(matrix, i);
      SCIP_Real lb = SCIPvarGetLbGlobal(var);
      SCIP_Real ub = SCIPvarGetUbGlobal(var);
      builder.setColLb(i, lb);
      builder.setColUb(i, ub);
      builder.setColLbInf(i, SCIPisInfinity(scip, -lb));
      builder.setColUbInf(i, SCIPisInfinity(scip, ub));

      builder.setColIntegral(i, SCIPvarIsIntegral(var));
      builder.setObj(i, SCIPvarGetObj(var));
   }

   /* set up rows */
   builder.setNumRows(nrows);
   for(int i = 0; i != nrows; ++i)
   {
      int* rowcols = SCIPmatrixGetRowIdxPtr(matrix, i);
      SCIP_Real* rowvals = SCIPmatrixGetRowValPtr(matrix, i);
      int rowlen = SCIPmatrixGetRowNNonzs(matrix, i);
      builder.addRowEntries(i, rowlen, rowcols, rowvals);

      SCIP_Real lhs = SCIPmatrixGetRowLhs(matrix, i);
      SCIP_Real rhs = SCIPmatrixGetRowRhs(matrix, i);
      builder.setRowLhs(i, lhs);
      builder.setRowRhs(i, rhs);
      builder.setRowLhsInf(i, SCIPisInfinity( scip, -lhs ));
      builder.setRowRhsInf(i, SCIPisInfinity( scip, rhs ));
   }

   return builder.build();
}

/*
 * Callback methods of presolver
 */

/** copy method for constraint handler plugins (called when SCIP copies plugins) */
static
SCIP_DECL_PRESOLCOPY(presolCopyMILP)
{  /*lint --e{715}*/
   SCIP_CALL( SCIPincludePresolMILP(scip) );

   return SCIP_OKAY;
}

/** destructor of presolver to free user data (called when SCIP is exiting) */
static
SCIP_DECL_PRESOLFREE(presolFreeMILP)
{  /*lint --e{715}*/
   SCIP_PRESOLDATA* data = SCIPpresolGetData(presol);
   assert(data != NULL);

   SCIPpresolSetData(presol, NULL);
   SCIPfreeBlockMemory(scip, &data);
   return SCIP_OKAY;
}

/** initialization method of presolver (called after problem was transformed) */
static
SCIP_DECL_PRESOLINIT(presolInitMILP)
{  /*lint --e{715}*/
   SCIP_PRESOLDATA* data = SCIPpresolGetData(presol);
   assert(data != NULL);

   data->lastncols = -1;
   data->lastnrows = -1;

   /* determine the number of threads to initialize the tbb scheduler with.
    * The tbb default is to use the number of available hardware threads */
   int numthreads = tbb::task_scheduler_init::default_num_threads();

   /* if data->threads has a lower number than that and is not set to automatic,
    * then use that lower number of threads */
   if( data->threads != 0 && data->threads < numthreads )
      numthreads = data->threads;

   new (&data->schedulerinit) tbb::task_scheduler_init(numthreads);

   return SCIP_OKAY;
}

/** deinitialization method of presolver (called before transformed problem is freed) */
static
SCIP_DECL_PRESOLEXIT(presolExitMILP)
{  /*lint --e{715}*/
   SCIP_PRESOLDATA* data = SCIPpresolGetData(presol);
   assert(data != NULL);

   /* deinitilize tbb scheduler by calling destructor of initialization object */
   data->schedulerinit.~task_scheduler_init();

   return SCIP_OKAY;
}


/** execution method of presolver */
static
SCIP_DECL_PRESOLEXEC(presolExecMILP)
{  /*lint --e{715}*/
   SCIP_Bool initialized;
   SCIP_Bool complete;
   SCIP_MATRIX* matrix;
   SCIP_PRESOLDATA* data;

   *result = SCIP_DIDNOTRUN;

   /* at the moment the presolve library does some weak dual reductions within the core */
   if( !SCIPallowWeakDualReds(scip) )
      return SCIP_OKAY;

   data = SCIPpresolGetData(presol);

   int nvars = SCIPgetNVars(scip);
   int nconss = SCIPgetNConss(scip);

   /* run only if the problem size reduced by some amount since the last call or if it is the first call */
   if( data->lastncols != -1 && data->lastnrows != -1 &&
       nvars > data->lastncols * 0.85 &&
       nconss > data->lastnrows * 0.85 )
      return SCIP_OKAY;

   SCIP_CALL( SCIPmatrixCreate(scip, &matrix, TRUE, &initialized, &complete) );

   /* we only work on pure MIPs */
   if( !initialized || !complete )
   {
      data->lastncols = 0;
      data->lastnrows = 0;

      if( initialized )
         SCIPmatrixFree(scip, &matrix);

      return SCIP_OKAY;
   }

   Problem<SCIP_Real> problem = buildProblem(scip, matrix);
   Presolve<SCIP_Real> presolve;

   /* important so that SCIP does not throw an error, e.g. when an integer variable is substituted
    * into a knapsack constraint */
   presolve.getPresolveOptions().substitutebinarieswithints = false;

   /* currently these changes cannot be communicated to SCIP correctly since a constraint needs
    * to be modified in the cases where slackvariables are removed from constraints but for the
    * presolve library those look like normal substitution on the postsolve stack */
   presolve.getPresolveOptions().removeslackvars = false;

   /* communicate the SCIP parameter to the presolve libary */
   presolve.getPresolveOptions().maxfillinpersubstitution = data->maxfillinpersubstitution;

   /* set up the presolvers that shall participate */
   using uptr = std::unique_ptr<PresolveMethod<SCIP_Real>>;

   presolve.addPresolveMethod( uptr( new CoefficientStrengthening<SCIP_Real>() ) );
   presolve.addPresolveMethod( uptr( new SimpleProbing<SCIP_Real>() ) );
   presolve.addPresolveMethod( uptr( new ConstraintPropagation<SCIP_Real>() ) );
   presolve.addPresolveMethod( uptr( new ImplIntDetection<SCIP_Real>() ) );
   presolve.addPresolveMethod( uptr( new FixContinuous<SCIP_Real>() ) );

   if( data->enableparallelrows )
      presolve.addPresolveMethod( uptr( new ParallelRowDetection<SCIP_Real>() ) );

   presolve.addPresolveMethod( uptr( new SimpleSubstitution<SCIP_Real>() ) );
   presolve.addPresolveMethod( uptr( new SimplifyInequalities<SCIP_Real>() ) );
   presolve.addPresolveMethod( uptr( new SingletonCols<SCIP_Real>() ) );
   presolve.addPresolveMethod( uptr( new DualFix<SCIP_Real>() ) );

   if( data->enablemultiaggr )
      presolve.addPresolveMethod( uptr( new Substitution<SCIP_Real>() ) );

   if( data->enableprobing )
      presolve.addPresolveMethod( uptr( new Probing<SCIP_Real>() ) );

   if( data->enablesparsify )
      presolve.addPresolveMethod( uptr( new Sparsify<SCIP_Real>() ) );

   if( data->enabledualinfer )
      presolve.addPresolveMethod( uptr( new DualInfer<SCIP_Real>() ) );

   /* domincated columns, stuffing, and if possible in SCIP the future, parallel columns
    * are strong dual reductions */
   if( SCIPallowStrongDualReds(scip) )
   {
      presolve.addPresolveMethod( uptr( new SingletonStuffing<SCIP_Real>() ) );

      if( data->enabledomcol )
         presolve.addPresolveMethod( uptr( new DominatedCols<SCIP_Real>() ) );

      /* todo: parallel cols cannot be handled by SCIP currently
       * addPresolveMethod( uptr( new ParallelColDetection<SCIP_Real>() ) ); */
   }

   /* set tolerances */
   presolve.setEpsilon(SCIPepsilon(scip));
   presolve.setFeasTol(SCIPfeastol(scip));

   /* adjust output settings of presolve libary */
#ifdef SCIP_PRESOLLIB_ENABLE_OUTPUT
   problem.setName(SCIPgetProbName(scip));
#else
   presolve.setVerbosityLevel(VerbosityLevel::QUIET);
#endif

   /* call the presolving */
   SCIPverbMessage(scip, SCIP_VERBLEVEL_HIGH, NULL,
               "   (%.1fs) running MILP presolver\n", SCIPgetSolvingTime(scip));
   int oldnnz = problem.getConstraintMatrix().getNnz();
   PresolveResult<SCIP_Real> res = presolve.apply(problem);
   data->lastncols = problem.getNCols();
   data->lastnrows = problem.getNRows();

   /* evaluate the result */
   switch(res.status)
   {
      case PresolveStatus::INFEASIBLE:
         *result = SCIP_CUTOFF;
         SCIPverbMessage(scip, SCIP_VERBLEVEL_HIGH, NULL,
               "   (%.1fs) MILP presolver detected infeasibility\n",
               SCIPgetSolvingTime(scip));
         SCIPmatrixFree(scip, &matrix);
         return SCIP_OKAY;
      case PresolveStatus::UNBND_OR_INFEAS:
      case PresolveStatus::UNBOUNDED:
         *result = SCIP_UNBOUNDED;
         SCIPverbMessage(scip, SCIP_VERBLEVEL_HIGH, NULL,
               "   (%.1fs) MILP presolver detected unboundedness\n",
               SCIPgetSolvingTime(scip));
         SCIPmatrixFree(scip, &matrix);
         return SCIP_OKAY;
      case PresolveStatus::UNCHANGED:
         *result = SCIP_DIDNOTFIND;
         data->lastncols = 0;
         data->lastnrows = 0;
         SCIPverbMessage(scip, SCIP_VERBLEVEL_HIGH, NULL,
               "   (%.1fs) MILP presolver found nothing\n",
               SCIPgetSolvingTime(scip));
         SCIPmatrixFree(scip, &matrix);
         return SCIP_OKAY;
      case PresolveStatus::REDUCED:
         data->lastncols = problem.getNCols();
         data->lastnrows = problem.getNRows();
         *result = SCIP_SUCCESS;
   }

   /* result indicated success, now populate the changes into the SCIP structures */
   std::vector<SCIP_VAR*> tmpvars;
   std::vector<SCIP_Real> tmpvals;

   /* if the number of nonzeros decreased by a sufficient factor, rather create all constraints from scratch */
   int newnnz = problem.getConstraintMatrix().getNnz();

   if( newnnz <= data->modifyconsfac * oldnnz )
   {
      int oldnrows = SCIPmatrixGetNRows(matrix);
      int newnrows = problem.getNRows();

      /* capture constraints that are still present in the problem after presolve */
      for( int i = 0; i < newnrows; ++i )
      {
         SCIP_CONS* c = SCIPmatrixGetCons(matrix, res.postsolve.origrow_mapping[i]);
         SCIP_CALL( SCIPcaptureCons(scip, c) );
      }

      /* delete all constraints */
      *ndelconss += oldnrows;
      *naddconss += newnrows;

      for( int i = 0; i < oldnrows; ++i )
      {
         SCIP_CALL( SCIPdelCons(scip, SCIPmatrixGetCons(matrix, i)) );
      }

      /* now loop over rows of presolved problem and create them as new linear constraints,
       * then release the old constraint after its name was passed to the new constraint */
      const Vec<RowFlags>& rflags = problem.getRowFlags();
      const auto& consmatrix = problem.getConstraintMatrix();
      for( int i = 0; i < newnrows; ++i )
      {
         auto rowvec = consmatrix.getRowCoefficients(i);
         const int* rowcols = rowvec.getIndices();
         /* SCIPcreateConsBasicLinear() requires a non const pointer */
         SCIP_Real* rowvals = const_cast<SCIP_Real*>(rowvec.getValues());
         int rowlen = rowvec.getLength();

         /* retrieve SCIP compatible left and right hand sides */
         SCIP_Real lhs = rflags[i].test(RowFlag::LHS_INF) ? - SCIPinfinity(scip) : consmatrix.getLeftHandSides()[i];
         SCIP_Real rhs = rflags[i].test(RowFlag::RHS_INF) ? SCIPinfinity(scip) : consmatrix.getRightHandSides()[i];

         /* create variable array matching the value array */
         tmpvars.clear();
         tmpvars.reserve(rowlen);
         for( int j = 0; j < rowlen; ++j )
            tmpvars.push_back(SCIPmatrixGetVar(matrix, res.postsolve.origcol_mapping[rowcols[j]]));

         /* create and add new constraint with name of old constraint */
         SCIP_CONS* oldcons = SCIPmatrixGetCons(matrix, res.postsolve.origrow_mapping[i]);
         SCIP_CONS* cons;
         SCIP_CALL( SCIPcreateConsBasicLinear(scip, &cons, SCIPconsGetName(oldcons), rowlen, tmpvars.data(), rowvals, lhs, rhs) );
         SCIP_CALL( SCIPaddCons(scip, cons) );

         /* release old and new constraint */
         SCIP_CALL( SCIPreleaseCons(scip, &oldcons) );
         SCIP_CALL( SCIPreleaseCons(scip, &cons) );
      }
   }

   /* loop over res.postsolve and add all fixed variables and aggregations to scip */
   for( std::size_t i = 0; i != res.postsolve.types.size(); ++i )
   {
      ReductionType type = res.postsolve.types[i];
      int first = res.postsolve.start[i];
      int last = res.postsolve.start[i + 1];

      switch( type )
      {
      case ReductionType::FIXED_COL:
      {
         SCIP_Bool infeas;
         SCIP_Bool fixed;
         int col = res.postsolve.indices[first];

         SCIP_VAR* colvar = SCIPmatrixGetVar(matrix, col);

         SCIP_Real value = res.postsolve.values[first];

         SCIP_CALL( SCIPfixVar(scip, colvar, value, &infeas, &fixed) );
         *nfixedvars += 1;

         assert(!infeas);
         assert(fixed);
         break;
      }
      case ReductionType::SUBSTITUTED_COL:
      {
         int col = res.postsolve.indices[first];
         SCIP_Real side = res.postsolve.values[first];
         SCIP_Real colCoef = 0.0;
         tmpvars.clear();
         tmpvals.clear();
         tmpvars.reserve(last - first - 1);
         tmpvals.reserve(last - first - 1);
         for( int j = first + 1; j < last; ++j )
         {
            if( res.postsolve.indices[j] == col )
            {
               colCoef = res.postsolve.values[j];
               break;
            }
         }

         assert(colCoef != 0.0);
         SCIP_VAR* aggrvar = SCIPmatrixGetVar(matrix, col);
         while( SCIPvarGetStatus(aggrvar) == SCIP_VARSTATUS_AGGREGATED )
         {
            SCIP_Real scalar = SCIPvarGetAggrScalar(aggrvar);
            SCIP_Real constant = SCIPvarGetAggrConstant(aggrvar);
            aggrvar = SCIPvarGetAggrVar(aggrvar);

            side -= colCoef * constant;
            colCoef *= scalar;
         }

         assert(SCIPvarGetStatus(aggrvar) != SCIP_VARSTATUS_MULTAGGR);

         for( int j = first + 1; j < last; ++j )
         {
            if( res.postsolve.indices[j] == col )
               continue;

            tmpvars.push_back(SCIPmatrixGetVar(matrix, res.postsolve.indices[j]));
            tmpvals.push_back(- res.postsolve.values[j] / colCoef);
         }

         SCIP_Bool infeas;
         SCIP_Bool aggregated;
         SCIP_CALL( SCIPmultiaggregateVar(scip, aggrvar, tmpvars.size(),
            tmpvars.data(), tmpvals.data(), side / colCoef, &infeas, &aggregated) );

         if( aggregated )
            *naggrvars += 1;

         if( infeas )
         {
            *result = SCIP_CUTOFF;
            break;
         }

         break;
      }
      default:
      case ReductionType::PARALLEL_COL:
         return SCIP_INVALIDRESULT;
      }
   }

   /* tighten bounds of variables that are still present after presolving */
   if( *result != SCIP_CUTOFF )
   {
      VariableDomains<SCIP_Real>& varDomains = problem.getVariableDomains();
      for( int i = 0; i != problem.getNCols(); ++i )
      {
         SCIP_VAR* var = SCIPmatrixGetVar(matrix, res.postsolve.origcol_mapping[i]);
         if( !varDomains.flags[i].test(ColFlag::LB_INF) )
         {
            SCIP_Bool infeas;
            SCIP_Bool tightened;
            SCIP_CALL( SCIPtightenVarLb(scip, var, varDomains.lower_bounds[i], TRUE, &infeas, &tightened) );

            if( tightened )
               *nchgbds += 1;

            if( infeas )
            {
               *result = SCIP_CUTOFF;
               break;
            }
         }

         if( !varDomains.flags[i].test(ColFlag::UB_INF) )
         {
            SCIP_Bool infeas;
            SCIP_Bool tightened;
            SCIP_CALL( SCIPtightenVarUb(scip, var, varDomains.upper_bounds[i], TRUE, &infeas, &tightened) );

            if( tightened )
               *nchgbds += 1;

            if( infeas )
            {
               *result = SCIP_CUTOFF;
               break;
            }
         }
      }
   }

   /* finish with a final verb message and return */
   SCIPverbMessage(scip, SCIP_VERBLEVEL_HIGH, NULL,
      "   (%.1fs) MILP presolver (%d rounds): %d aggregations, %d fixings, %d bound changes\n",
      SCIPgetSolvingTime(scip), presolve.getStatistics().nrounds, *naggrvars, *nfixedvars, *nchgbds);

   /* free the matrix */
   assert(initialized);
   SCIPmatrixFree(scip, &matrix);

   return SCIP_OKAY;
}


/*
 * presolver specific interface methods
 */

/** creates the xyz presolver and includes it in SCIP */
SCIP_RETCODE SCIPincludePresolMILP(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_PRESOLDATA* presoldata;
   SCIP_PRESOL* presol;

   /* create MILP presolver data */
   presoldata = NULL;
   SCIP_CALL( SCIPallocBlockMemory(scip, &presoldata) );

   presol = NULL;

   /* include presolver */
   SCIP_CALL( SCIPincludePresolBasic(scip, &presol, PRESOL_NAME, PRESOL_DESC, PRESOL_PRIORITY, PRESOL_MAXROUNDS, PRESOL_TIMING,
         presolExecMILP,
         presoldata) );

   assert(presol != NULL);

   /* set non fundamental callbacks via setter functions */
   SCIP_CALL( SCIPsetPresolCopy(scip, presol, presolCopyMILP) );
   SCIP_CALL( SCIPsetPresolFree(scip, presol, presolFreeMILP) );
   SCIP_CALL( SCIPsetPresolInit(scip, presol, presolInitMILP) );
   SCIP_CALL( SCIPsetPresolExit(scip, presol, presolExitMILP) );

   /* add MILP presolver parameters */
   SCIP_CALL( SCIPaddIntParam(scip,
         "presolving/" PRESOL_NAME "/threads",
         "maximum number of threads presolving may use (0: automatic)",
         &presoldata->threads, FALSE, DEFAULT_THREADS, 0, INT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(scip,
         "presolving/" PRESOL_NAME "/maxfillinpersubstitution",
         "maximal possible fillin for substitutions to be considered",
         &presoldata->maxfillinpersubstitution, FALSE, DEFAULT_MAXFILLINPERSUBST, INT_MIN, INT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip,
         "presolving/" PRESOL_NAME "/modifyconsfac",
         "modify SCIP constraints when the number of nonzeros is at most this factor times the number of nonzeros before presolving",
         &presoldata->modifyconsfac, FALSE, DEFAULT_MODIFYCONSFAC, 0.0, 1.0, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip,
         "presolving/" PRESOL_NAME "/enableparallelrows",
         "should the parallel rows presolver be enabled within the presolve library?",
         &presoldata->enableparallelrows, TRUE, DEFAULT_ENABLEPARALLELROWS, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip,
         "presolving/" PRESOL_NAME "/enabledomcol",
         "should the dominated column presolver be enabled within the presolve library?",
         &presoldata->enabledomcol, TRUE, DEFAULT_ENABLEDOMCOL, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip,
         "presolving/" PRESOL_NAME "/enabledualinfer",
         "should the dualinfer presolver be enabled within the presolve library?",
         &presoldata->enabledualinfer, TRUE, DEFAULT_ENABLEDUALINFER, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip,
         "presolving/" PRESOL_NAME "/enablemultiaggr",
         "should the multi-aggregation presolver be enabled within the presolve library?",
         &presoldata->enablemultiaggr, TRUE, DEFAULT_ENABLEMULTIAGGR, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip,
         "presolving/" PRESOL_NAME "/enableprobing",
         "should the probing presolver be enabled within the presolve library?",
         &presoldata->enableprobing, TRUE, DEFAULT_ENABLEPROBING, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip,
         "presolving/" PRESOL_NAME "/enablesparsify",
         "should the sparsify presolver be enabled within the presolve library?",
         &presoldata->enablesparsify, TRUE, DEFAULT_ENABLESPARSIFY, NULL, NULL) );

   return SCIP_OKAY;
}

#endif
