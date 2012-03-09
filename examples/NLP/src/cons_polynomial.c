/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2012 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   cons_polynomial.h
 * @brief  constraint handler for polynomial constraints
 * @author Lignfeng Niu, 
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>

#include "scip/pub_lp.h"
#include "scip/struct_var.h"
#include "cons_polynomial.h"

#include "scip/pub_var.h"

#include "solve_nlp.h"

/* constraint handler properties */
#define CONSHDLR_NAME          "polynomial"
#define CONSHDLR_DESC          "polynomial constraints"
#define CONSHDLR_SEPAPRIORITY         0 /**< priority of the constraint handler for separation */
#define CONSHDLR_ENFOPRIORITY  -4000000 /**< priority of the constraint handler for constraint enforcing */
#define CONSHDLR_CHECKPRIORITY -4000000 /**< priority of the constraint handler for checking feasibility */
#define CONSHDLR_SEPAFREQ            -1 /**< frequency for separating cuts; zero means to separate only in the root node */
#define CONSHDLR_PROPFREQ            -1 /**< frequency for propagating domains; zero means only preprocessing propagation */
#define CONSHDLR_EAGERFREQ          100 /**< frequency for using all instead of only the useful constraints in separation,
                                              *   propagation and enforcement, -1 for no eager evaluations, 0 for first only */
#define CONSHDLR_MAXPREROUNDS        -1 /**< maximal number of presolving rounds the constraint handler participates in (-1: no limit) */
#define CONSHDLR_DELAYSEPA        FALSE /**< should separation method be delayed, if other separators found cuts? */
#define CONSHDLR_DELAYPROP        FALSE /**< should propagation method be delayed, if other propagators found reductions? */
#define CONSHDLR_DELAYPRESOL      FALSE /**< should presolving method be delayed, if other presolvers found reductions? */
#define CONSHDLR_NEEDSCONS         TRUE /**< should the constraint handler be skipped, if no constraints are available? */

/* TODO: (optional) enable linear constraint upgrading */
#if 0
#include "scip/cons_linear.h"

#define LINCONSUPGD_PRIORITY          0
#endif




/*
 * Data structures
 */

/* TODO: fill in the necessary constraint data */


struct SCIP_ConsData
{
   Polynomial* polynomial;
   SCIP_Real lhs;
   SCIP_Real rhs;
};

/** constraint handler data */
struct SCIP_ConshdlrData
{
};

/*
 * Static data used in this file
 */

/* ??? */
static SCIP_VAR ** scip_vars = NULL;
/* ??? */
static int scip_nvars, scip_nbinvars, scip_nintvars, scip_nimplvars, scip_ncontvars, scip_begin_of_contvars;


/*
 * Local methods
 */

/* put your local methods here, and declare them static */

/* get the 0-based index of a var */
static 
int indexOfVar( SCIP* scip, SCIP_VAR* var )
{
   int i;
   if( scip_vars == NULL )
   { /* Load the vars if not loaded */
      SCIP_CALL( SCIPgetVarsData(scip, &scip_vars, &scip_nvars, &scip_nbinvars, &scip_nintvars, &scip_nimplvars, &scip_ncontvars) );
   }

      
  assert( scip_vars != NULL );

  scip_begin_of_contvars = scip_nbinvars + scip_nintvars + scip_nimplvars;
  assert( scip_begin_of_contvars + scip_ncontvars == scip_nvars );

  printf("the var need to be check the index is %s %p \n", SCIPvarGetName( var ), var);

  printf("scip_begin_of_contvars is %d, and scip_nvars is %d, scip_ncontvars is %d \n", scip_begin_of_contvars, scip_nvars, scip_ncontvars);
  for( i = scip_begin_of_contvars; i < scip_nvars; ++i )
  {
    /* compare the variable by address based on the fact that one variable has only one copy in memory */
    printf("the current one in scip_vars is %s %p \n", SCIPvarGetName( scip_vars[i] ), scip_vars[i]);

    if( scip_vars[i] == var ) return i; 
  }
  printf("not find var %f", var->obj);
  abort();
  /**  assert( 0 );  terminate the program if run here, since no index for var found */
  return -1;
}

/* evaluate the value of a monomial at a point given in SCIP_SOL or the LP soltuion (sol = NULL) */
static
SCIP_Real evaluateMonomial( SCIP * scip , SCIP_SOL * sol , Monomial * monomial )
{
   int i;
   SCIP_Real val = monomial->coefficient;
   printf("The monomial's coefficient is %g \n", monomial->coefficient);

   /* ??? Why can val be 0 ? */

   for( i = 0; i < monomial->nvars && val != 0; ++i)
   {
      val = val * pow( SCIPgetSolVal(scip, sol, monomial->vars[i]), monomial->power[i] );
      printf("The monomial var <%s>=%g, power=%g, val=%g", SCIPvarGetName(monomial->vars[i]),SCIPgetSolVal(scip, sol, monomial->vars[i]), monomial->power[i],val);
   }

   return val;
}

/* evaluate the value of a polynomial at a point given in SCIP_SOL or the LP soltuion (sol = NULL) */
static
SCIP_Real evaluatePolynomial(SCIP * scip, SCIP_SOL * sol, Polynomial * polynomial )
{
   int i;
   SCIP_Real val = 0;
   for( i = 0; i < polynomial->nMonomials ; ++i )
   {
      val = val + evaluateMonomial(scip, sol, polynomial->monomials[i]);
   }
   return val;
}

/** given the constriant structure and a solution or the LP soltuion (sol = NULL), and evaluate this constriant
 *  value at this special point */
static
SCIP_Real evaluateConstraint(SCIP * scip, SCIP_SOL * sol, SCIP_CONSDATA * consdata)
{
        SCIP_Real val;
   val = evaluatePolynomial(scip, sol, consdata->polynomial);
   return val;
}

/** checks polynomial constraint for feasibility of given solution or the LP soltuion (sol = NULL), returns TRUE iff constraint is violated */
static
SCIP_RETCODE checkConstraintViolated(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint to check */
   SCIP_SOL*             sol,                /**< solution to check, NULL for current solution */
   SCIP_Bool*            violated            /**< pointer to store whether the constraint is violated */
   )
{
   SCIP_Real val;
   SCIP_CONSDATA* consdata;

   assert(violated != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   SCIPdebugMessage("checking polynomial constraint <%s> for feasibility of solution %p \n",
      SCIPconsGetName(cons), sol);

   *violated = FALSE;
   val = evaluateConstraint(scip, sol, consdata);
   
   printf("The current polynomial is:\n");
   {
         int mono;
         int v;
         /* print left hand side for ranged rows */
         if( !SCIPisInfinity(scip, -consdata->lhs)
            && !SCIPisInfinity(scip, consdata->rhs)
            && !SCIPisEQ(scip, consdata->lhs, consdata->rhs) )
            SCIPinfoMessage(scip, NULL, "%g <= ", consdata->lhs);

         /* print the polynomial */
         for( mono  = 0; mono < consdata->polynomial->nMonomials; ++mono)
         {
            assert(consdata->polynomial->monomials[mono] != NULL );

            /** print the coefficient for the current monomial */
            SCIPinfoMessage(scip, NULL , "%+g", consdata->polynomial->monomials[mono]->coefficient);

            /** print the variables and powers for the current monomial */
            for( v=0; v < consdata->polynomial->monomials[mono]->nvars; ++v)
            {
               SCIPinfoMessage(scip, NULL , "pow(<%s>,%f)", SCIPvarGetName(consdata->polynomial->monomials[mono]->vars[v]), consdata->polynomial->monomials[mono]->power[v]);
            }
         }
   }


   printf("val=%g, lhs=%g, rhs = %g\n", val, consdata->rhs, consdata->lhs);

   /*   if(val > consdata->rhs || val < consdata->lhs)*/
   if( SCIPisFeasGT( scip, val, consdata->rhs) || SCIPisFeasLT( scip, val, consdata->lhs))
   {
      *violated = TRUE;
      return SCIP_OKAY;
   }
   return SCIP_OKAY;
}

/** there are some local methods to create the data structure step by step */

/** create a monomial */
SCIP_RETCODE monomialCreate(
   SCIP*                 scip,                /**< SCIP data structure */
   Monomial**            monomials,           /**< pointer to store constraint data */
   int                   nvars,  
   SCIP_VAR**            vars, 
   SCIP_Real*            power,
   SCIP_Real             coefficient 
   )
{
   assert( monomials != NULL );

   SCIP_CALL( SCIPallocBlockMemory(scip, monomials) );

   (*monomials)->nvars = nvars;
   (*monomials)->coefficient = coefficient;
   SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, & ( (*monomials)->vars ), vars , nvars) );
   SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, & ( (*monomials)->power), power, nvars ) );

   (*monomials)->nuses = 1;
   {
      int i;
      puts("Hiho");
      for(i = 0; i < (*monomials)->nvars; i++)
         printf("%d %g %g\n", i, (*monomials)->power[i],power[i]);
   }
   return SCIP_OKAY;   
}

static
SCIP_RETCODE monomialFree(
   SCIP*         scip, 
   Monomial**    monomials
   )
{
   assert( monomials != NULL);
   assert( *monomials != NULL );
   SCIPfreeBlockMemoryArray(scip, &( (*monomials)->vars ), (*monomials)->nvars);
   SCIPfreeBlockMemoryArray(scip, &( (*monomials)->power) ,(*monomials)->nvars);
   SCIPfreeBlockMemory(scip, monomials);
   *monomials = NULL;

   return SCIP_OKAY;
}

static
SCIP_Bool is_monomial_valid(
   Monomial *monomial
   )
{
   int i;
   
   assert( monomial != NULL);
   assert( monomial->nuses >= 1);
   assert( monomial->nvars >= 1);

   for(i=0; i< monomial->nvars; ++i)
   {
      if( monomial->power[i] == 0 )
      {
         printf("monomial should not have power 0!");
         return FALSE;
      }
   }
   return TRUE;
}

static
void captureMonomial(
   SCIP*      scip,
        Monomial*  monomial
   )
{
   assert(monomial != NULL);
   monomial->nuses ++;
}

void releaseMonomial(
   SCIP*        scip, 
   Monomial**   monomials
   )
{
   assert(  monomials );
   assert( *monomials );
   (*monomials)->nuses --;
   if( (*monomials)->nuses == 0 )
      monomialFree(scip, monomials );
   *monomials = NULL;
}
    

/** create a polynomial */
SCIP_RETCODE polynomialCreate(
   SCIP*          scip, 
   Polynomial**   polynomials,
   int            nMonomials, 
   Monomial**     monomials)
{
   int i;
   assert( polynomials != NULL );
   assert( monomials != NULL );

   SCIP_CALL( SCIPallocBlockMemory(scip, polynomials) );
   (*polynomials)->nMonomials = nMonomials;
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &( (*polynomials)->monomials), nMonomials ) );
   for(i = 0; i < nMonomials; ++i)
   {
      assert( monomials[i] != NULL );
      (*polynomials)->monomials[i] = monomials[i];
      captureMonomial(scip, (*polynomials)->monomials[i] );
   }
   (*polynomials)->nuses = 1;
   return SCIP_OKAY;
}

static
SCIP_Bool is_polynomial_valid(
   Polynomial *polynomial
   )
{
   int i;
   
   assert( polynomial != NULL);
   assert( polynomial->nuses >= 1);
   assert( polynomial->nMonomials >= 1);
   
   for(i=0; i<polynomial->nMonomials; ++i)
   {
      if(is_monomial_valid(polynomial->monomials[i]) != TRUE)
      {
         printf("Invalid monomial member makes Polynomial fail!");
         return FALSE;
      }
   }
   
   return TRUE;
}


static
void polynomialFree(
   SCIP*         scip, 
   Polynomial**  polynomials)
{
   int i;
   assert( polynomials );
   assert( *polynomials);

   for(i = 0; i < (*polynomials)->nMonomials; ++i)
   {
        releaseMonomial(scip,&( (*polynomials)->monomials[i]) );
   }
   assert((*polynomials)->monomials != NULL);
   
   SCIPfreeBlockMemoryArray(scip, &(*polynomials)->monomials, (*polynomials)->nMonomials );
   SCIPfreeBlockMemory(scip, polynomials);
   *polynomials = NULL;
}

static
void capturePolynomial(
   SCIP*             scip, 
   Polynomial*       polynomial
   )
{
   assert(polynomial != NULL);
   polynomial->nuses++;
}

void releasePolynomial(
   SCIP*          scip, 
   Polynomial**   polynomials)
{
   assert( polynomials );
   assert( *polynomials );
   (*polynomials)->nuses--;
   if( (*polynomials)->nuses == 0 )
      polynomialFree(scip, polynomials );
   *polynomials = NULL;
}

/** creates polynomial constraint data */
static
SCIP_RETCODE consdataCreate(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSDATA**       consdata,           /**< pointer to store constraint data */
   Polynomial*           polynomial,
   SCIP_Real             lhs,
   SCIP_Real             rhs
   )
{
   int i;
   int j;
   assert(consdata != NULL);
   /**   assert(*consdata != NULL); */
   assert(polynomial != NULL);

   SCIP_CALL( SCIPallocBlockMemory(scip, consdata) );
   (*consdata)->lhs = lhs; 
   (*consdata)->rhs = rhs;
   (*consdata)->polynomial = polynomial;
    capturePolynomial(scip, (*consdata)->polynomial);

   for( i=0; i<(*consdata)->polynomial->nMonomials; ++i)
   {
      for(j=0; j<(*consdata)->polynomial->monomials[i]->nvars; ++j)
      {
         /* get transformed variable */
         printf("name %s %p \n", SCIPvarGetName( (*consdata)->polynomial->monomials[i]->vars[j] ), (*consdata)->polynomial->monomials[i]->vars[j]);
      }
   }
  
   if( SCIPisTransformed(scip) )
   {
      printf("after construct the consdata, we need transform it\n");
      for( i=0; i<(*consdata)->polynomial->nMonomials; ++i)
      {
         for(j=0; j<(*consdata)->polynomial->monomials[i]->nvars; ++j)
         {
            /* get transformed variable */
            SCIP_CALL(
               SCIPgetTransformedVar( scip, (*consdata)->polynomial->monomials[i]->vars[j], &((*consdata)->polynomial->monomials[i]->vars[j])) );
            printf("name %s %p \n", SCIPvarGetName( (*consdata)->polynomial->monomials[i]->vars[j] ), (*consdata)->polynomial->monomials[i]->vars[j]); 
         }
      }
   }

   return SCIP_OKAY;
}

static
SCIP_Bool is_consdata_valid(
   SCIP_CONSDATA *consdata
   )
{
   assert( consdata != NULL );
   
   if(is_polynomial_valid(consdata->polynomial) != TRUE)
   {
      printf("Invalid polynomial makes consdata fail!");
      return FALSE;
   }

   return TRUE;
}
/** prints polynomial constraint to file stream */
static
void consdataPrint(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSDATA*        consdata,           /**< polynomial constraint data */
   FILE*                 file                /**< output file (or NULL for standard output) */
   )
{
   int mono;
   int v;

   assert(consdata != NULL);
   assert(consdata->polynomial != NULL );
   
   /* print left hand side for ranged rows */
   if( !SCIPisInfinity(scip, -consdata->lhs)
      && !SCIPisInfinity(scip, consdata->rhs)
      && !SCIPisEQ(scip, consdata->lhs, consdata->rhs) )
      SCIPinfoMessage(scip, file, "%g <= ", consdata->lhs);

   /* print the polynomial */
   for( mono  = 0; mono < consdata->polynomial->nMonomials; ++mono)
   {
       assert(consdata->polynomial->monomials[mono] != NULL );

      /** print the coefficient for the current monomial */
      SCIPinfoMessage(scip, file, "%+g", consdata->polynomial->monomials[mono]->coefficient);
      
      /** print the variables and powers for the current monomial */
      for( v=0; v < consdata->polynomial->monomials[mono]->nvars; ++v)
      {
         SCIPinfoMessage(scip, file, "pow(<%s>,%f)", SCIPvarGetName(consdata->polynomial->monomials[mono]->vars[v]), consdata->polynomial->monomials[mono]->power[v]);         
      }
   }
   
   /* print right hand side */
   if( SCIPisEQ(scip, consdata->lhs, consdata->rhs) )
      SCIPinfoMessage(scip, file, "== %g\n", consdata->rhs);
   else if( !SCIPisInfinity(scip, consdata->rhs) )
      SCIPinfoMessage(scip, file, "<= %g\n", consdata->rhs);
   else if( !SCIPisInfinity(scip, -consdata->lhs) )
      SCIPinfoMessage(scip, file, ">= %g\n", consdata->lhs);
   else
      SCIPinfoMessage(scip, file, " [free]\n");
}

static 
SCIP_RETCODE ipoptSolve( SCIP* scip, SCIP_CONS** conss, int nconss, SCIP_SOL* sol, SCIP_Real * contVarsVals )
{
   int n_lprows = 0;
   int n_lpcols = 0;
   int i;
   int j;
   int k;
   int ind;
   int nvars;
   SCIP_VAR** vars = NULL;
   SCIP_ROW** rows = NULL;
   SCIP_COL** cols = NULL;
   SCIP_Real val;
   SCIP_CONSDATA* consdata;
   NLP nlp;

   int        tem_int;
   SCIP_Real  tem_real;
   SCIP_Real* tem_real_array = NULL;
   SCIP_VAR*  tem_var        = NULL;
   SCIP_VAR** tem_var_array  = NULL;

   /* prepare for the NLP structure for Ipopt in four steps */
   
   /* prepare some parameters for alloc memoery and fill in the simple data in NLP structure */

   /* since we already called SCIP_CALL( SCIPgetVarsData(scip, &scip_vars, &scip_nvars, &scip_nbinvars, &scip_nimplvars, &scip_ncontvars) ); 
    *  in the improvedByiIpopt, so in scip_vars should already instore the pointer for vars, we only need to check here */
   assert(scip_vars != NULL );

   printf("gather the general information for NLP. \n");
   //   nlp.nvars = scip_ncontvars;
   nlp.nvars = scip_nvars;
   nlp.nbinvars = scip_nbinvars;
   nlp.nintvars = scip_nintvars;
   nlp.nimplvars = scip_nimplvars;
   nlp.ncontvars = scip_ncontvars;
   nlp.nactivevars = SCIPgetNVars(scip);
   nlp.nnonactivevars = SCIPgetNFixedVars(scip);

   assert(nlp.nvars == nlp.nbinvars + nlp.nintvars + nlp.nimplvars + nlp.ncontvars);
   assert(nlp.nvars == nlp.nactivevars + nlp.nnonactivevars);

   SCIP_CALL( SCIPgetLPRowsData (scip, &rows, &n_lprows) );

   nlp.m_LP  = n_lprows + nlp.nnonactivevars;
   nlp.m_NLP = nconss;
   //   nlp.mcons = nconss+n_lprows;
   nlp.mcons = nlp.m_LP + nlp.m_NLP;

   printf("There are %d variables, %d constraints, %d are linear constraits, %d are nonlinear constriants.\n", nlp.nvars, nlp.mcons, nlp.m_LP, nlp.m_NLP);
   
   SCIP_CALL( SCIPallocBufferArray( scip, &(nlp.c),  nlp.nvars ) );
   SCIP_CALL( SCIPallocBufferArray( scip, &(nlp.x),  nlp.nvars ) );
   SCIP_CALL( SCIPallocBufferArray( scip, &(nlp.lb), nlp.nvars ) );
   SCIP_CALL( SCIPallocBufferArray( scip, &(nlp.ub), nlp.nvars ) );

   SCIP_CALL( SCIPallocBufferArray( scip, &(nlp.nnonz),  nlp.m_LP ) );
   SCIP_CALL( SCIPallocBufferArray( scip, &(nlp.jCols),  nlp.m_LP ) );
   SCIP_CALL( SCIPallocBufferArray( scip, &(nlp.values), nlp.m_LP ) );

   SCIP_CALL( SCIPallocBufferArray( scip, &(nlp.polynomials), nlp.m_NLP ) );

   SCIP_CALL( SCIPallocBufferArray( scip, &(nlp.rhs), nlp.mcons ) ); 
   SCIP_CALL( SCIPallocBufferArray( scip, &(nlp.lhs), nlp.mcons ) );


   /** prepare the data for linear objective function and variables */
   for( j = 0; j < nlp.nvars; ++j )
   {
      nlp.lb[j] = - SCIP_DEFAULT_INFINITY; 
      nlp.ub[j] =   SCIP_DEFAULT_INFINITY; /* its value is 10^20 */
   }

   /** since when scip call ipopt, all the integer include binary aready fixed to integer values,
       we fix their lb and ub the integer values. Notice these variables might also occur in the
       NLP part
   */
   for(j = 0; j< (nlp.nbinvars + nlp.nintvars) ; ++j)
   {
      nlp.c[j]  = SCIPvarGetObj(scip_vars[j] );
      nlp.x[j]  = SCIPgetSolVal( scip, sol, scip_vars[j] );
      nlp.lb[j] = nlp.x[j];
      nlp.ub[j] = nlp.x[j];
   }
   
   /**@todo now, we do not fix the implicite integral part, since we also get the LP part from scip,
      the soultion should be integer for there vairables, if not, either the presolving of scip is
      wrong, or the solution given by ipopt is not exact for these variables
   */
   for(j =(nlp.nbinvars + nlp.nintvars) ; j<nlp.nvars; ++j)
   {
      nlp.c[j]  = SCIPvarGetObj(scip_vars[j] );
      nlp.x[j]  = SCIPgetSolVal( scip, sol, scip_vars[j] );
      nlp.lb[j] = SCIPvarGetLbLocal( scip_vars[j] );
      nlp.ub[j] = SCIPvarGetUbLocal( scip_vars[j] );      
   }
      
   printf("start to get the LP part in ipoptSol");
   SCIP_CALL( SCIPgetLPColsData(scip, &cols, &n_lpcols) );
   for( i = 0; i < n_lprows; ++i )
   {
       val = SCIProwGetConstant( rows[i] );

       nlp.rhs[i] = SCIProwGetRhs( rows[i] ) - val;
       nlp.lhs[i] = SCIProwGetLhs( rows[i] ) - val;
       nlp.nnonz[i] = SCIProwGetNLPNonz( rows[i] );

       SCIP_CALL( SCIPallocBufferArray( scip, &( nlp.jCols[i] ),nlp.nnonz[i] ) );
        
       SCIP_CALL( SCIPduplicateBufferArray( scip, &( nlp.values[i] ), SCIProwGetVals( rows[i] ), nlp.nnonz[i] ) );

       for( j = 0; j < nlp.nnonz[i] ; ++j )
       {
       /* get the index of var */
       printf("name %s %p \n", SCIPvarGetName( SCIPcolGetVar( cols[j] )), SCIPcolGetVar( cols[j] ));
       ind = indexOfVar( scip, SCIPcolGetVar( cols[j] ) );
       //ind = ind - scip_begin_of_contvars;
       /* set the coefficient in row sparse mode */
       nlp.jCols[i][j] = ind;
       }
   }

   printf("\n the LP part finished.\n");
   
   /** reconstruct the missing part because of presolving */
   assert(i == n_lprows);
   
   nlp.nactivevars = SCIPgetNVars(scip);
   nlp.nnonactivevars = SCIPgetNFixedVars(scip);
   nlp.nfixed = 0;
   nlp.naggr = 0;
   nlp.nmultaggr = 0;
   nlp.nnegation = 0;

   if(nlp.nnonactivevars == 0)
   {
      vars = SCIPgetFixedVars(scip);

      for( k=0; k < nlp.nnonactivevars ; ++k)
      {
         /** for each fixed variable, we add a linear constriant to it */

         switch( SCIPvarGetStatus( vars[k] ))
         {
         case SCIP_VARSTATUS_FIXED:     /** 3 variable is fixed to specific value in the transformed problem */
            nlp.rhs[i] =  SCIPvarGetLbLocal(vars[k]);
            nlp.lhs[i] =  SCIPvarGetLbLocal(vars[k]);
            nlp.nnonz[i] = 1;
            SCIP_CALL( SCIPallocBufferArray( scip, &( nlp.jCols[i] ),nlp.nnonz[i] ) );
            SCIP_CALL( SCIPallocBufferArray( scip, &( nlp.values[i]),nlp.nnonz[i] ) );

            nlp.values[i][0] = 1.0;
            printf("name %s %p \n", SCIPvarGetName( vars[k] ), vars[k] );
            ind = indexOfVar( scip, vars[k] );
            //ind = ind - scip_begin_of_contvars;
            nlp.jCols[i][0] = ind;
            nlp.nfixed += 1;
            i += 1;

            break;
         case SCIP_VARSTATUS_AGGREGATED: /** 4  variable is aggregated to x = a*y + c in the transformed problem */
            tem_var = SCIPvarGetAggrVar( vars[k] );
            tem_real = SCIPvarGetAggrConstant( vars[k] );
            nlp.rhs[i] = tem_real;
            nlp.lhs[i] = tem_real;
            nlp.nnonz[i] = 2;
            SCIP_CALL( SCIPallocBufferArray( scip, &( nlp.jCols[i] ),nlp.nnonz[i] ) );
            SCIP_CALL( SCIPallocBufferArray( scip, &( nlp.values[i]),nlp.nnonz[i] ) );

            tem_real = SCIPvarGetAggrScalar( vars[k] );

            nlp.values[i][0] = 1.0;
            printf("name %s %p \n", SCIPvarGetName( vars[k] ), vars[k] );
            ind = indexOfVar( scip, vars[k] );
            //ind = ind - scip_begin_of_contvars;
            nlp.jCols[i][0] = ind;

            nlp.values[i][1] = -tem_real;
            printf("name %s %p \n", SCIPvarGetName( tem_var ), tem_var );
            ind = indexOfVar( scip, tem_var);
            //ind = ind - scip_begin_of_contvars;
            nlp.jCols[i][1] = ind;
            nlp.naggr += 1;
            i += 1;

            break;
         case SCIP_VARSTATUS_MULTAGGR:   /** 5  variable is aggregated to x = a_1*y_1 + ... + a_k*y_k + c */
            tem_real = SCIPvarGetMultaggrConstant(vars[k]);
            nlp.lhs[i] = -tem_real;
            nlp.rhs[i] = -tem_real;

            tem_int = SCIPvarGetMultaggrNVars( vars[k] );
            nlp.nnonz[i] = tem_int + 1;

            SCIP_CALL( SCIPallocBufferArray( scip, &( nlp.jCols[i] ),nlp.nnonz[i] ) );
            SCIP_CALL( SCIPallocBufferArray( scip, &( nlp.values[i]),nlp.nnonz[i] ) );
            tem_var_array = SCIPvarGetMultaggrVars( vars[k] );
            tem_real_array =  SCIPvarGetMultaggrScalars( vars[k] );

            for( j = 0; j < tem_int; ++j )
            {
               nlp.values[i][j] = tem_real_array[j];
               printf("name %s %p \n", SCIPvarGetName( tem_var_array[j] ), tem_var_array[j] );
               ind = indexOfVar( scip, tem_var_array[j] );
               //ind = ind - scip_begin_of_contvars;
               nlp.jCols[i][j] = ind;
            }

            nlp.values[i][j] = -1.0;
            printf("name %s %p \n", SCIPvarGetName( vars[k] ), vars[k] );
            ind = indexOfVar( scip, vars[k] );
            //ind = ind - scip_begin_of_contvars;
            nlp.jCols[i][j] = ind;
            nlp.nmultaggr += 1;
            i += 1;

            break;
         case SCIP_VARSTATUS_NEGATED:  /** 6  variable is the negation of an original or transformed variable */
            /** gets the negation variable x of a negated variable x' = offset - x */
            tem_real = SCIPvarGetNegationConstant( vars[k] );
            nlp.lhs[i] = tem_real;
            nlp.rhs[i] = tem_real;

            tem_var = SCIPvarGetNegationVar( vars[k]);

            nlp.nnonz[i] = 2;
            SCIP_CALL( SCIPallocBufferArray( scip, &( nlp.jCols[i] ),nlp.nnonz[i] ) );
            SCIP_CALL( SCIPallocBufferArray( scip, &( nlp.values[i]),nlp.nnonz[i] ) );

            nlp.values[i][0] = 1.0;
            printf("name %s %p \n", SCIPvarGetName( vars[k] ), vars[k] );
            ind = indexOfVar( scip, vars[k] );
            //ind = ind - scip_begin_of_contvars;
            nlp.jCols[i][0] = ind;

            nlp.values[i][1] = 1.0;
            printf("name %s %p \n", SCIPvarGetName( tem_var ), tem_var );
            ind = indexOfVar( scip, tem_var );
            // ind = ind - scip_begin_of_contvars;
            nlp.jCols[i][1] = ind;

            i += 1;
            break;
         default:

            printf("var %p is not a var deliminate in the presolving, its status is %d, please check code.\n", vars[j],SCIPvarGetStatus( vars[k] ));
            exit(0);
         }
      }
   }

   /** prepare data for polynomial nonlinear constraints */

   for( i = 0; i < nlp.m_NLP ; ++i )  /* for each polynomial in constraints */
   {
       consdata = SCIPconsGetData(conss[i]);
       assert(consdata != NULL);
       assert( is_consdata_valid(consdata) == TRUE);
       consdataPrint(scip,consdata,NULL);
   
       SCIP_CALL( SCIPallocBuffer(scip, &(nlp.polynomials[i]) ));
       nlp.polynomials[i]->nMonomials = consdata->polynomial->nMonomials;
       SCIP_CALL( SCIPallocBufferArray( scip, &(nlp.polynomials[i]->monomials), nlp.polynomials[i]->nMonomials ) );

       for( j = 0; j < nlp.polynomials[i]->nMonomials; ++j ) /* for each monomial */
       {
          SCIP_CALL( SCIPallocBuffer( scip, & (nlp.polynomials[i]->monomials[j]) ));
          nvars =  consdata->polynomial->monomials[j]->nvars;
          nlp.polynomials[i]->monomials[j]->nvars = nvars;
          nlp.polynomials[i]->monomials[j]->coefficient = consdata->polynomial->monomials[j]->coefficient;

          SCIP_CALL( SCIPallocBufferArray(scip, &(nlp.polynomials[i]->monomials[j]->indicies ), nvars ));
          SCIP_CALL( SCIPallocBufferArray(scip, &(nlp.polynomials[i]->monomials[j]->power ), nvars ));

          printf("now we print each variable's name and address!\n");
          for( k = 0; k < nvars; ++k ) /* for each variable */
          {
             printf("name %s %p \n", SCIPvarGetName( consdata->polynomial->monomials[j]->vars[k] ), consdata->polynomial->monomials[j]->vars[k]);
          }
          printf("print this monomial part is finished! \n");

          
          for( k = 0; k < nvars; ++k ) /* for each variable */
          {
             printf("name %s %p \n", SCIPvarGetName( consdata->polynomial->monomials[j]->vars[k] ), consdata->polynomial->monomials[j]->vars[k]);

             printf("we get the index and power of this variable\n");
             // nlp.polynomials[i]->monomials[j]->indicies[k] = indexOfVar( scip, consdata->polynomial->monomials[j]->vars[k] ) - scip_begin_of_contvars;
             nlp.polynomials[i]->monomials[j]->indicies[k] = indexOfVar( scip, consdata->polynomial->monomials[j]->vars[k] );
             nlp.polynomials[i]->monomials[j]->power[k] = consdata->polynomial->monomials[j]->power[k];
             printf("we finish the job of getting the index and power for this variable\n");
             /* need to count the variables in the polynomial here */
             /* buffer[monomials[j]->indicies[k]] = 1; */
          }
          printf("we finish the job of getting data for the current monomial\n");
       }

       nlp.lhs[i + nlp.m_LP] = consdata->lhs;
       nlp.rhs[i + nlp.m_LP] = consdata->rhs;
      
   }
   /* until now we already fill all the information for nlp */


   printf("There are %d variables, %d constraints, %d are linear constraits, %d are nonlinear constriants.\n", nlp.nvars, nlp.mcons, nlp.m_LP, nlp.m_NLP);


   /** print the nlp struct here, let us see whether we get all the information, and whethere they are correct*/

   PrintNLP( scip, &nlp );
   
   /* call ipopt to get a improved solution */
   Callipopt( scip, &nlp );

   PrintNLP( scip, &nlp);

   /** pass the solution back to improvedByIpopt */
   for(i=0; i<nlp.nvars ; ++i)
   {
      contVarsVals[i] = nlp.x[i];
   }
   /** free the memory for nlp here */

   SCIPfreeBufferArray( scip, &(nlp.c) );
   SCIPfreeBufferArray( scip, &(nlp.x) );
   SCIPfreeBufferArray( scip, &(nlp.ub) );
   SCIPfreeBufferArray( scip, &(nlp.lb) );
   for(i=0; i<nlp.m_LP; ++i)
   {
       SCIPfreeBufferArray( scip, &( nlp.jCols[i]  ) );
       SCIPfreeBufferArray( scip, &( nlp.values[i] ) );
   }
   SCIPfreeBufferArray( scip, &(nlp.values));
   SCIPfreeBufferArray( scip, &(nlp.jCols));
   SCIPfreeBufferArray( scip, &(nlp.nnonz) );
   for(i=0; i<nlp.m_NLP; ++i)
   {
      for(j=0; j<nlp.polynomials[i]->nMonomials; ++j)
      {
         SCIPfreeBufferArray( scip, &(nlp.polynomials[i]->monomials[j]->indicies ));
         SCIPfreeBufferArray( scip, &(nlp.polynomials[i]->monomials[j]->power ));
         SCIPfreeBuffer( scip, &(nlp.polynomials[i]->monomials[j]));
      }
     
      SCIPfreeBufferArray( scip, &(nlp.polynomials[i]->monomials));
      SCIPfreeBuffer( scip, &(nlp.polynomials[i]));
   }
   SCIPfreeBufferArray(scip, &(nlp.polynomials));
   SCIPfreeBufferArray( scip, &(nlp.lhs) );
   SCIPfreeBufferArray( scip, &(nlp.rhs) );
   
   /* SCIPfreeBuffer( scip, nlp ) ; */

   /*SCIPfreeBufferArray( &buffer );*/
   return SCIP_OKAY;
}

/** in this function we get a SCIP_SOL or a LP solution as a input (notice we already check this point do not satisfiy 
 * the whole polynomial constriants, otherwise we can not arrive there), call ipopt to comupte a improved solution and 
 * store it in the SCIP solution pool. I am not sure now, it is the mechanism I thought SCIP works.*/
static
SCIP_RETCODE improveSolByIpopt(
   SCIP *      scip,
   SCIP_CONSHDLR * conshdlr,
   SCIP_CONS ** conss,
   int nconss,
   SCIP_SOL* sol
   )
{
   int i;
   SCIP_SOL * improvedSol;
   SCIP_Bool  stored;
   SCIP_Real * contVarsVals;

   
   /* create primal solution */
   printf("call SCIPcreateSol in improveSolByIpopt to get a new sol to be improved \n");

   SCIP_CALL( SCIPcreateSol(scip, &improvedSol, NULL) );

   /* since IPOPT only improved the solution of the continuous part, we get the values of binary and integer variables from the current node */
   printf("call SCIPgetVarsData to get data for the linear part and store it in the global vars");

   SCIP_CALL( SCIPgetVarsData(scip, &scip_vars, &scip_nvars, &scip_nbinvars, &scip_nintvars, &scip_nimplvars, &scip_ncontvars) );

   
   
   for( i = 0; i< scip_begin_of_contvars ; ++i )
   {
      SCIP_CALL( SCIPsetSolVal(scip, improvedSol ,scip_vars[i],  SCIPgetSolVal(scip, sol , scip_vars[i]))); /* ???? */
   } 
   
     /** improvedSol = sol; */
   SCIPallocBufferArray(scip,&contVarsVals,scip_ncontvars);
   /* now we call ipopt to improve the continuous part */
   printf(" now we call ipopt to improve the continuous part ");

   ipoptSolve(scip, conss, nconss, sol, contVarsVals);

   
   printf("get the improved part and pass it back to SCIP\n");
   for( i = scip_begin_of_contvars; i < scip_nvars; ++i )
   {
      /* set the solution value of the variable */
      /* SCIP_CALL( SCIPsetSolVal(scip, improvedSol ,scip_vars[i], nlp->x[i- scip_begin_of_contvars]) );*/
      SCIP_CALL( SCIPsetSolVal(scip, improvedSol ,scip_vars[i], contVarsVals[i-scip_begin_of_contvars]) );
      printf("x<%d> should be %g, is %g \n", i-scip_begin_of_contvars, contVarsVals[i-scip_begin_of_contvars],
         SCIPgetSolVal(scip, improvedSol, scip_vars[i])); 
   }

   /** or */
   /* SCIP_CALL( SCIPsetSolVals(scip, improvedSol, nvars, vars, values); */   
        /** add solution to the solution storage and free solution afterwards */
   SCIPfreeBufferArray(scip,&contVarsVals);
   printf("old solution is %p with transobj %g and origobj %g\n adding solution %p with transobj %g and origobj %g\n",
     sol, SCIPgetSolTransObj(scip, sol), SCIPgetSolOrigObj(scip, sol),
      improvedSol, SCIPgetSolTransObj(scip, improvedSol), SCIPgetSolOrigObj(scip, improvedSol));
   SCIP_CALL( SCIPaddSolFree(scip, &improvedSol, &stored) );
   return SCIP_OKAY;
}


/*
 * Linear constraint upgrading
 */

#ifdef LINCONSUPGD_PRIORITY
/** tries to upgrade a linear constraint into a polynomial constraint */
static
SCIP_DECL_LINCONSUPGD(linconsUpgdPolynomial)
{  /*lint --e{715}*/
   SCIP_Bool upgrade;

   assert(upgdcons != NULL);
   
   /* check, if linear constraint can be upgraded to polynomial constraint */
   upgrade = FALSE;
   /* TODO: put the constraint's properties here, in terms of the statistics given by nposbin, nnegbin, ... */

   if( upgrade )
   {
      SCIPdebugMessage("upgrading constraint <%s> to polynomial constraint\n", SCIPconsGetName(cons));
      
      /* create the bin Polynomial constraint (an automatically upgraded constraint is always unmodifiable) */
      assert(!SCIPconsIsModifiable(cons));
      SCIP_CALL( SCIPcreateConsPolynomial(scip, upgdcons, SCIPconsGetName(cons), nvars, vars, vals, lhs, rhs,
            SCIPconsIsInitial(cons), SCIPconsIsSeparated(cons), SCIPconsIsEnforced(cons), 
            SCIPconsIsChecked(cons), SCIPconsIsPropagated(cons), SCIPconsIsLocal(cons),
            SCIPconsIsLocal(cons), SCIPconsIsModifiable(cons), 
            SCIPconsIsDynamic(cons), SCIPconsIsRemovable(cons), SCIPconsIsStickingAtNode(cons)) );
   }

   return SCIP_OKAY;
}
#endif



/*
 * Callback methods of constraint handler
 */

/* TODO: Implement all necessary constraint handler methods. The methods with #if 0 ... #else #define ... are optional */

/** destructor of constraint handler to free constraint handler data (called when SCIP is exiting) */
#if 0
static
SCIP_DECL_CONSFREE(consFreePolynomial)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of polynomial constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consFreePolynomial NULL
#endif


/** initialization method of constraint handler (called after problem was transformed) */
#if 1
static
SCIP_DECL_CONSINIT(consInitPolynomial)
{  /*lint --e{715}*/
   assert(scip != null);

   if( SCIPfindConshdlr(scip,"linear") != NULL )
   {
#ifdef LINCONSUPGD_PRIORITY
      /* include the linear constraint upgrade in the linear constraint handler */
      SCIP_CALL( SCIPincludeLinconsUpgrade(scip, linconsUpgdPolynomial, LINCONSUPGD_PRIORITY, CONSHDLR_NAME) );
#endif
   }
   return SCIP_OKAY;
}
#else
#define consInitPolynomial NULL
#endif


/** deinitialization method of constraint handler (called before transformed problem is freed) */
#if 0
static
SCIP_DECL_CONSEXIT(consExitPolynomial)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of polynomial constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consExitPolynomial NULL
#endif


/** presolving initialization method of constraint handler (called when presolving is about to begin) */
#if 0
static
SCIP_DECL_CONSINITPRE(consInitprePolynomial)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of polynomial constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consInitprePolynomial NULL
#endif


/** presolving deinitialization method of constraint handler (called after presolving has been finished) */
#if 0
static
SCIP_DECL_CONSEXITPRE(consExitprePolynomial)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of polynomial constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consExitprePolynomial NULL
#endif


/** solving process initialization method of constraint handler (called when branch and bound process is about to begin) */
#if 0
static
SCIP_DECL_CONSINITSOL(consInitsolPolynomial)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of polynomial constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consInitsolPolynomial NULL
#endif


/** solving process deinitialization method of constraint handler (called before branch and bound process data is freed) */
#if 0
static
SCIP_DECL_CONSEXITSOL(consExitsolPolynomial)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of polynomial constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consExitsolPolynomial NULL
#endif


/** frees specific constraint data */
#if 0
static
SCIP_DECL_CONSDELETE(consDeletePolynomial)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of polynomial constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consDeletePolynomial NULL
#endif


/** transforms constraint data into data belonging to the transformed problem */ 
#if 1
static
SCIP_DECL_CONSTRANS(consTransPolynomial)
{
   /** SCIP_CONSHDLRDATA* conshdlrdata; */
   SCIP_CONSDATA*     sourcedata;
   SCIP_CONSDATA*     targetdata;

   /*debugMessage("Trans method of linear constraints\n");*/

   assert(conshdlr != NULL);
   assert(SCIPgetStage(scip) == SCIP_STAGE_TRANSFORMING);
   assert(sourcecons != NULL);
   assert(targetcons != NULL);

   sourcedata = SCIPconsGetData(sourcecons);
   assert(sourcedata != NULL);
   assert(sourcedata->polynomial->monomials != NULL);

   /* create polynomial constraint data for target constraint */
   SCIP_CALL( consdataCreate(scip, &targetdata, sourcedata->polynomial, sourcedata->lhs, sourcedata->rhs) );

   /* create target constraint */
   SCIP_CALL( SCIPcreateCons(scip, targetcons, SCIPconsGetName(sourcecons), conshdlr, targetdata,
         SCIPconsIsInitial(sourcecons), SCIPconsIsSeparated(sourcecons), SCIPconsIsEnforced(sourcecons),
         SCIPconsIsChecked(sourcecons), SCIPconsIsPropagated(sourcecons),
         SCIPconsIsLocal(sourcecons), SCIPconsIsModifiable(sourcecons),
         SCIPconsIsDynamic(sourcecons), SCIPconsIsRemovable(sourcecons), SCIPconsIsStickingAtNode(sourcecons)) );

   printf("after transforming:\n");
   consdataPrint(scip, targetdata, NULL);
   
   return SCIP_OKAY;

   


}
#else
#define consTransPolynomial NULL
#endif

/** LP initialization method of constraint handler */
#if 0
static
SCIP_DECL_CONSINITLP(consInitlpPolynomial)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of polynomial constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consInitlpPolynomial NULL
#endif


/** separation method of constraint handler for LP solutions */
#if 0
static
SCIP_DECL_CONSSEPALP(consSepalpPolynomial)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of polynomial constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consSepalpPolynomial NULL   
#endif


/** separation method of constraint handler for arbitrary primal solutions */
#if 0
static
SCIP_DECL_CONSSEPASOL(consSepasolPolynomial)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of polynomial constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consSepasolPolynomial NULL
#endif


/** constraint enforcing method of constraint handler for LP solutions */
static
SCIP_DECL_CONSENFOLP(consEnfolpPolynomial)
{  
   SCIP_Bool violated;
   int i;
   *result = SCIP_FEASIBLE;

   /* conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL); */

   /* search for violated polynomial constraints */
   for( i = 0; i < nconss; i++ )
   {
      SCIP_CALL( checkConstraintViolated(scip, conss[i], NULL, &violated) );
      if( violated )
      {
         *result = SCIP_INFEASIBLE;
         printf("this call improveSolByIpopt from consEnfolpPolynomial");

         SCIP_CALL(improveSolByIpopt(scip, conshdlr, conss, nconss, NULL));

         /* if all integers are fixed, no branching is possible and we
          * have to prune the node
          */
         if( SCIPgetNPseudoBranchCands(scip) == 0 )
            *result = SCIP_CUTOFF;

            return SCIP_OKAY;      
      }
   } 

   
   return SCIP_OKAY;
}


/** constraint enforcing method of constraint handler for pseudo solutions */
static
SCIP_DECL_CONSENFOPS(consEnfopsPolynomial)
{ 
   /*lint --e{715}*/
   SCIP_Bool violated;
   int i;
   SCIP_Bool cutoff;

   cutoff = FALSE;
   SCIP_CALL( SCIPconstructLP(scip, &cutoff) );
   if( cutoff == TRUE )
   {
      *result = SCIP_CUTOFF;
      return SCIP_OKAY;
   }
   else
   {
      for( i = 0; i < nconss; i++ )
      {
         SCIP_CALL( checkConstraintViolated(scip, conss[i], NULL, &violated) );
         if( violated )
         {
            *result = SCIP_INFEASIBLE;
            printf("this call improveSolByIpopt from consEnfopsPolynomial");

            SCIP_CALL(improveSolByIpopt(scip, conshdlr, conss, nconss, NULL));

            /* if all integers are fixed, no branching is possible and we
             * have to prune the node
             */
            if( SCIPgetNPseudoBranchCands(scip) == 0 )
               *result = SCIP_CUTOFF;
            return SCIP_OKAY;
         }
      } 
      *result = SCIP_FEASIBLE;
      return SCIP_OKAY;
   }
}


/** feasibility check method of constraint handler for integral solutions */
static
SCIP_DECL_CONSCHECK(consCheckPolynomial)
{ 
   /*lint --e{715}*/
   SCIP_Bool violated;
   int i;

   if(SCIPsolGetOrigin(sol) == SCIP_SOLORIGIN_ORIGINAL)
   {
      *result = SCIP_FEASIBLE;
      return SCIP_OKAY;
   }

   printf("checking polynomial constraint for solution %p with transobj %g and origobj %g\n",
      sol, SCIPgetSolTransObj(scip, sol), SCIPgetSolOrigObj(scip, sol));
   assert(sol != NULL);
   
   for( i = 0; i < nconss; i++ )
   {
      SCIP_CALL( checkConstraintViolated(scip, conss[i], sol, &violated) );
      if( violated )
      {
         *result = SCIP_INFEASIBLE;
         printf("this callimproveSolByIpopt is from consCheckPolynomial\n");

         if( SCIPgetStage(scip) <= SCIP_STAGE_SOLVING )
         {
            SCIP_CALL(improveSolByIpopt(scip, conshdlr, conss, nconss, sol));
         }
         return SCIP_OKAY;
      }
   } 
   *result = SCIP_FEASIBLE;

   return SCIP_OKAY;
   
   /**   for( i = 0; i < nconss; i++ )
   {
      SCIP_CALL( checkConstraintViolated(scip, conss[i], sol, &violated) );
      if( violated )
      {
         *result = SCIP_INFEASIBLE;
         return SCIP_OKAY;
      }
   } 
   *result = SCIP_FEASIBLE;

   return SCIP_OKAY; */
}


/** domain propagation method of constraint handler */
#if 0
static
SCIP_DECL_CONSPROP(consPropPolynomial)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of polynomial constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consPropPolynomial NULL
#endif


/** presolving method of constraint handler */
#if 0
static
SCIP_DECL_CONSPRESOL(consPresolPolynomial)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of polynomial constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consPresolPolynomial NULL
#endif


/** propagation conflict resolving method of constraint handler */
#if 0
static
SCIP_DECL_CONSRESPROP(consRespropPolynomial)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of polynomial constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consRespropPolynomial NULL
#endif


/** variable rounding lock method of constraint handler */
static
SCIP_DECL_CONSLOCK(consLockPolynomial)
{  /*lint --e{715}*/
   SCIP_CONSDATA* consdata;
   int i;
   int j;

   printf("Getting consdata\n");
   
   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);
   consdataPrint(scip, consdata, NULL);   /**    consdataPrint(scip, consdata, NULL); */
   for( i = 0; i < consdata->polynomial->nMonomials; ++i)
   {
      for( j = 0; j < consdata->polynomial->monomials[i]->nvars; ++j)
      {
         SCIP_CALL( SCIPaddVarLocks(scip, consdata->polynomial->monomials[i]->vars[j], nlockspos + nlocksneg, nlockspos + nlocksneg) );
      }
   }
   return SCIP_OKAY;
}


/** constraint activation notification method of constraint handler */
#if 0
static
SCIP_DECL_CONSACTIVE(consActivePolynomial)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of polynomial constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consActivePolynomial NULL
#endif


/** constraint deactivation notification method of constraint handler */
#if 0
static
SCIP_DECL_CONSDEACTIVE(consDeactivePolynomial)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of polynomial constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consDeactivePolynomial NULL
#endif


/** constraint enabling notification method of constraint handler */
#if 0
static
SCIP_DECL_CONSENABLE(consEnablePolynomial)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of polynomial constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consEnablePolynomial NULL
#endif


/** constraint disabling notification method of constraint handler */
#if 0
static
SCIP_DECL_CONSDISABLE(consDisablePolynomial)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of polynomial constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consDisablePolynomial NULL
#endif

/** variable deletion method of constraint handler */
#define consDelVarsPolynomial NULL

/** constraint display method of constraint handler */
#if 1
static
SCIP_DECL_CONSPRINT(consPrintPolynomial)
{  /*lint --e{715}*/
   consdataPrint(scip, SCIPconsGetData(cons), file);
   
   return SCIP_OKAY;
}
#else
#define consPrintPolynomial NULL
#endif


/** constraint copying method of constraint handler */
#if 0
static
SCIP_DECL_CONSCOPY(consCopyPolynomial)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of polynomial constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consCopyPolynomial NULL
#endif


/** constraint parsing method of constraint handler */
#if 0
static
SCIP_DECL_CONSPARSE(consParsePolynomial)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of polynomial constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consParsePolynomial NULL
#endif


/*
 * constraint specific interface methods
 */

/** creates the handler for polynomial constraints and includes it in SCIP */
SCIP_DECL_INCLUDEPLUGIN(SCIPincludeConshdlrPolynomial)
{
   SCIP_CONSHDLRDATA* conshdlrdata;

   /* create polynomial constraint handler data */
   conshdlrdata = NULL;
   /* TODO: (optional) create constraint handler specific data here */

   /* include constraint handler */
   SCIP_CALL( SCIPincludeConshdlr(scip, CONSHDLR_NAME, CONSHDLR_DESC,
         CONSHDLR_SEPAPRIORITY, CONSHDLR_ENFOPRIORITY, CONSHDLR_CHECKPRIORITY,
         CONSHDLR_SEPAFREQ, CONSHDLR_PROPFREQ, CONSHDLR_EAGERFREQ, CONSHDLR_MAXPREROUNDS, 
         CONSHDLR_DELAYSEPA, CONSHDLR_DELAYPROP, CONSHDLR_DELAYPRESOL, CONSHDLR_NEEDSCONS,
         SCIPincludeConshdlrPolynomial,
         consFreePolynomial, consInitPolynomial, consExitPolynomial, 
         consInitprePolynomial, consExitprePolynomial, consInitsolPolynomial, consExitsolPolynomial,
         consDeletePolynomial, consTransPolynomial, consInitlpPolynomial,
         consSepalpPolynomial, consSepasolPolynomial, consEnfolpPolynomial, consEnfopsPolynomial, consCheckPolynomial, 
         consPropPolynomial, consPresolPolynomial, consRespropPolynomial, consLockPolynomial,
         consActivePolynomial, consDeactivePolynomial, 
         consEnablePolynomial, consDisablePolynomial, consDelVarsPolynomial,
         consPrintPolynomial, consCopyPolynomial, consParsePolynomial,
         conshdlrdata) );

   /* add polynomial constraint handler parameters */
   /* TODO: (optional) add constraint handler specific parameters with SCIPaddTypeParam() here */

   return SCIP_OKAY;
}

/** creates and captures a polynomial constraint */
SCIP_RETCODE SCIPcreateConsPolynomial(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS**           cons,               /**< pointer to hold the created constraint */
   const char*           name,               /**< name of constraint */
   Polynomial*           polynomial,
   SCIP_Real             lhs,                /**< left hand side of constraint */
   SCIP_Real             rhs,                /**< right hand side of constraint */
   SCIP_Bool             initial,            /**< should the LP relaxation of constraint be in the initial LP? */
   SCIP_Bool             separate,           /**< should the constraint be separated during LP processing? */
   SCIP_Bool             enforce,            /**< should the constraint be enforced during node processing? */
   SCIP_Bool             check,              /**< should the constraint be checked for feasibility? */
   SCIP_Bool             propagate,          /**< should the constraint be propagated during node processing? */
   SCIP_Bool             local,              /**< is constraint only valid locally? */
   SCIP_Bool             modifiable,         /**< is constraint modifiable (subject to column generation)? */
   SCIP_Bool             dynamic,            /**< is constraint subject to aging? */
   SCIP_Bool             removable,          /**< should the relaxation be removed from the LP due to aging or cleanup? */
   SCIP_Bool             stickingatnode      /**< should the node always be kept at the node where it was added, even
                                              *   if it may be moved to a more global node? */
   )
{
   /* TODO: (optional) modify the definition of the SCIPcreateConsPolynomial() call, if you don't need all the information */
   SCIP_CONSHDLR* conshdlr;
   SCIP_CONSDATA* consdata;

   printf("Creating consdata\n");
   
   /* find the polynomial constraint handler */
   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);
   if( conshdlr == NULL )
   {
      SCIPerrorMessage("polynomial constraint handler not found\n");
      return SCIP_PLUGINNOTFOUND;
   }

   /* create and store constraint specific data here */

   SCIP_CALL( consdataCreate(scip, &consdata, polynomial, lhs, rhs) );

   assert(consdata != NULL);
   consdataPrint(scip, consdata, NULL);
   
   /* create constraint */
   SCIP_CALL( SCIPcreateCons(scip, cons, name, conshdlr, consdata, initial, separate, enforce, check, propagate, local, modifiable, dynamic, removable, stickingatnode) );

   assert(SCIPconsGetData(*cons) != NULL);

   assert(SCIPconsGetData(*cons) == consdata);
   assert(consdata->polynomial->monomials != NULL);
   
   printf("This is the constriants by calling SCIPconsGetData()\n");
   consdataPrint(scip, SCIPconsGetData(*cons), NULL);
   
   return SCIP_OKAY;
}
