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
 * @author Robert Lion Gottwald
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>

#include "scip/presol_milp.h"
#include "core/Presolve.hpp"
#include "core/ProblemBuilder.hpp"


#define PRESOL_NAME            "milp"
#define PRESOL_DESC            "MILP specific presolving routine"
#define PRESOL_PRIORITY        -9999999 /**< priority of the presolver (>= 0: before, < 0: after constraint handlers); combined with propagators */
#define PRESOL_MAXROUNDS             -1 /**< maximal number of presolving rounds the presolver participates in (-1: no limit) */
#define PRESOL_TIMING           SCIP_PRESOLTIMING_MEDIUM /* timing of the presolver (fast, medium, or exhaustive) */


/*
 * Data structures
 */

/* TODO: fill in the necessary presolver data */

/** presolver data */
struct SCIP_PresolData
{
   int lastncols;
   int lastnrows;
};


/*
 * Local methods
 */

/* put your local methods here, and declare them static */


/*
 * Callback methods of presolver
 */

/* TODO: Implement all necessary presolver methods. The methods with an #if 0 ... #else #define ... are optional */


/** copy method for constraint handler plugins (called when SCIP copies plugins) */
static
SCIP_DECL_PRESOLCOPY(presolCopyMILP)
{  /*lint --e{715}*/
   SCIP_CALL( SCIPincludePresolMILP(scip) );

   return SCIP_OKAY;
}

/** destructor of presolver to free user data (called when SCIP is exiting) */
#if 1
static
SCIP_DECL_PRESOLFREE(presolFreeMILP)
{  /*lint --e{715}*/
   SCIP_PRESOLDATA* data = SCIPpresolGetData(presol);
   assert(data != NULL);

   SCIPpresolSetData(presol, NULL);
   SCIPfreeBlockMemory(scip, &data);
   return SCIP_OKAY;
}
#else
#define presolFreeMILP NULL
#endif

/** initialization method of presolver (called after problem was transformed) */
#if 1
static
SCIP_DECL_PRESOLINIT(presolInitMILP)
{  /*lint --e{715}*/
   SCIP_PRESOLDATA* data = SCIPpresolGetData(presol);
   assert(data != NULL);

   data->lastncols = -1;
   data->lastnrows = -1;

   return SCIP_OKAY;
}
#else
#define presolInitMILP NULL
#endif


/** deinitialization method of presolver (called before transformed problem is freed) */
#if 0
static
SCIP_DECL_PRESOLEXIT(presolExitMILP)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of xyz presolver not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define presolExitMILP NULL
#endif


/** presolving initialization method of presolver (called when presolving is about to begin) */
#if 0
static
SCIP_DECL_PRESOLINITPRE(presolInitpreMILP)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of xyz presolver not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define presolInitpreMILP NULL
#endif


/** presolving deinitialization method of presolver (called after presolving has been finished) */
#if 0
static
SCIP_DECL_PRESOLEXITPRE(presolExitpreMILP)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of xyz presolver not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define presolExitpreMILP NULL
#endif

static
Problem<SCIP_Real>
buildProblem(SCIP* scip, SCIP_MATRIX* matrix)
{
   ProblemBuilder<SCIP_Real> builder;

   // build problem from matrix
   int nnz = SCIPmatrixGetNNonzs(matrix);
   int ncols = SCIPmatrixGetNColumns(matrix);
   int nrows = SCIPmatrixGetNRows(matrix);
   builder.reserve(nnz, nrows, ncols);
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

/** execution method of presolver */
static
SCIP_DECL_PRESOLEXEC(presolExecMILP)
{  /*lint --e{715}*/
   SCIP_Bool initialized;
   SCIP_Bool complete;
   SCIP_MATRIX* matrix;
   SCIP_PRESOLDATA* data;

   *result = SCIP_DIDNOTRUN;
   if( SCIPgetNRuns(scip) != 1 )
      return SCIP_OKAY;

   data = SCIPpresolGetData(presol);

   int nvars = SCIPgetNVars(scip);
   int nconss = SCIPgetNConss(scip);

   if( data->lastncols != -1 && data->lastnrows != -1 &&
       nvars > data->lastncols * 0.85 &&
       nconss > data->lastnrows * 0.85 )
      return SCIP_OKAY;

   SCIP_CALL( SCIPmatrixCreate(scip, &matrix, TRUE, &initialized, &complete) );

   if( !initialized || !complete )
   {
      data->lastncols = 0;
      data->lastnrows = 0;

      if( initialized )
         SCIPmatrixFree(scip, &matrix);

      return SCIP_OKAY;
   }

   /* we only work on pure MIPs */
   Problem<SCIP_Real> problem = buildProblem(scip, matrix);
   Presolve<SCIP_Real> presolve;

   using uptr = std::unique_ptr<PresolveMethod<SCIP_Real>>;

   presolve.addPresolveMethod( uptr( new CoefficientStrengthening<SCIP_Real>() ) );
   presolve.addPresolveMethod( uptr( new SimpleProbing<SCIP_Real>() ) );
   presolve.addPresolveMethod( uptr( new ConstraintPropagation<SCIP_Real>() ) );
   presolve.addPresolveMethod( uptr( new ImplIntDetection<SCIP_Real>() ) );
   presolve.addPresolveMethod( uptr( new FixContinuous<SCIP_Real>() ) );
   presolve.addPresolveMethod( uptr( new ParallelRowDetection<SCIP_Real>() ) );
   // todo: parallel cols cannot be handled by SCIP currently
   // addPresolveMethod( uptr( new ParallelColDetection<SCIP_Real>() ) );
   presolve.addPresolveMethod( uptr( new SimpleSubstitution<SCIP_Real>() ) );
   presolve.addPresolveMethod( uptr( new Substitution<SCIP_Real>() ) );
   presolve.addPresolveMethod( uptr( new Probing<SCIP_Real>() ) );
   // presolve.addPresolveMethod( uptr( new Sparsify<SCIP_Real>() ) );
   presolve.addPresolveMethod( uptr( new SimplifyInequalities<SCIP_Real>() ) );

   if( SCIPallowWeakDualReds(scip) )
   {
      presolve.addPresolveMethod( uptr( new SingletonCols<SCIP_Real>() ) );
      presolve.addPresolveMethod( uptr( new DualFix<SCIP_Real>() ) );
      presolve.addPresolveMethod( uptr( new DualInfer<SCIP_Real> ) );
   }

   if( SCIPallowStrongDualReds(scip) )
   {
      presolve.addPresolveMethod( uptr( new SingletonStuffing<SCIP_Real>() ) );
      presolve.addPresolveMethod( uptr( new DominatedCols<SCIP_Real>() ) );
   }

   presolve.setEpsilon(SCIPepsilon(scip));
   presolve.setFeasTol(SCIPfeastol(scip));
   presolve.setVerbosityLevel(VerbosityLevel::QUIET);

   PresolveResult<SCIP_Real> res = presolve.apply(problem);
   data->lastncols = problem.getNCols();
   data->lastnrows = problem.getNRows();

   switch(res.status)
   {
      case PresolveStatus::INFEASIBLE:
         *result = SCIP_CUTOFF;
         SCIPmatrixFree(scip, &matrix);
         return SCIP_OKAY;
      case PresolveStatus::UNBOUNDED:
         *result = SCIP_UNBOUNDED;
         SCIPmatrixFree(scip, &matrix);
         return SCIP_OKAY;
      case PresolveStatus::UNBND_OR_INFEAS:
         //todo
      case PresolveStatus::UNCHANGED:
         *result = SCIP_DIDNOTFIND;
         data->lastncols = 0;
         data->lastnrows = 0;
         puts("didnotfind");
         SCIPmatrixFree(scip, &matrix);
         return SCIP_OKAY;
      case PresolveStatus::REDUCED:
         data->lastncols = problem.getNCols();
         data->lastnrows = problem.getNRows();
         puts("success");
         *result = SCIP_SUCCESS;
   }

   std::vector<SCIP_VAR*> aggrvars;
   std::vector<SCIP_Real> aggrvals;

   // tighten bounds of variables that are still present
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

   // loop over res.postsolve and add all fixed variables and aggregations to scip
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
         aggrvars.clear();
         aggrvals.clear();
         aggrvars.reserve(last - first - 1);
         aggrvals.reserve(last - first - 1);
         for( int j = first + 1; j < last; ++j )
         {
            if( res.postsolve.indices[j] == col )
            {
               colCoef = res.postsolve.values[j];
               break;
            }
         }

         assert(colCoef != 0.0);

         for( int j = first + 1; j < last; ++j )
         {
            if( res.postsolve.indices[j] == col )
               continue;

            aggrvars.push_back(SCIPmatrixGetVar(matrix, res.postsolve.indices[j]));
            aggrvals.push_back(- res.postsolve.values[j] / colCoef);
         }

         SCIP_Bool infeas;
         SCIP_Bool aggregated;
         SCIP_CALL( SCIPmultiaggregateVar(scip, SCIPmatrixGetVar(matrix, col), aggrvars.size(),
            aggrvars.data(), aggrvals.data(),side / colCoef, &infeas, &aggregated) );

         if( aggregated )
         {
            *naggrvars += 1;
         }

         if( infeas )
         {
            *result = SCIP_CUTOFF;
            break;
         }

         break;
      }
      case ReductionType::PARALLEL_COL:
         assert(false);
      }
   }

   if( initialized )
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

   /* create xyz presolver data */
   presoldata = NULL;
   /* TODO: (optional) create presolver specific data here */
   SCIP_CALL( SCIPallocBlockMemory(scip, &presoldata) );

   presol = NULL;

   /* include presolver */
#if 0
   /* use SCIPincludePresol() if you want to set all callbacks explicitly and realize (by getting compiler errors) when
    * new callbacks are added in future SCIP versions
    */
   SCIP_CALL( SCIPincludePresol(scip, PRESOL_NAME, PRESOL_DESC, PRESOL_PRIORITY, PRESOL_MAXROUNDS, PRESOL_TIMING,
         presolCopyXyz, presolFreeXyz, presolInitXyz, presolExitXyz, presolInitpreXyz, presolExitpreXyz, presolExecXyz,
         presoldata) );
#else
   /* use SCIPincludePresolBasic() plus setter functions if you want to set callbacks one-by-one and your code should
    * compile independent of new callbacks being added in future SCIP versions
    */
   SCIP_CALL( SCIPincludePresolBasic(scip, &presol, PRESOL_NAME, PRESOL_DESC, PRESOL_PRIORITY, PRESOL_MAXROUNDS, PRESOL_TIMING,
         presolExecMILP,
         presoldata) );

   assert(presol != NULL);

   /* set non fundamental callbacks via setter functions */
   SCIP_CALL( SCIPsetPresolCopy(scip, presol, presolCopyMILP) );
   SCIP_CALL( SCIPsetPresolFree(scip, presol, presolFreeMILP) );
   SCIP_CALL( SCIPsetPresolInit(scip, presol, presolInitMILP) );
   SCIP_CALL( SCIPsetPresolExit(scip, presol, presolExitMILP) );
   SCIP_CALL( SCIPsetPresolInitpre(scip, presol, presolInitpreMILP) );
   SCIP_CALL( SCIPsetPresolExitpre(scip, presol, presolExitpreMILP) );
#endif

   /* add MILP presolver parameters */
   /* TODO: (optional) add presolver specific parameters with SCIPaddTypeParam() here */

   return SCIP_OKAY;
}
