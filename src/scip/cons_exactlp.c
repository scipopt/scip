/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2009 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#pragma ident "@(#) $Id: cons_exactlp.c,v 1.1.2.1 2009/07/13 12:48:48 bzfwolte Exp $"
// #define SCIP_DEBUG /*??????????????*/
// #define CONSEXLP_OUT /* only for debugging ???????????????? */

/**@file   cons_exactlp.c
 * @ingroup CONSHDLRS 
 * @brief  constraint handler for exactlp constraints
 * @author Tobias Achterberg
 * @author Kati Wolter
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>
#include <gmp.h> 

#include "scip/cons_exactlp.h"
#include "EGlib.h" 
#include "QSopt_ex.h" 


/* constraint handler properties */
#define CONSHDLR_NAME          "exactlp"
#define CONSHDLR_DESC          "LP relaxation of a MIP that is given by rational data"
#define CONSHDLR_SEPAPRIORITY    950000 /**< priority of the constraint handler for separation */
#define CONSHDLR_ENFOPRIORITY   -400000 /**< priority of the constraint handler for constraint enforcing */
#define CONSHDLR_CHECKPRIORITY  -400000 /**< priority of the constraint handler for checking feasibility */
#define CONSHDLR_SEPAFREQ             1 /**< frequency for separating cuts; zero means to separate only in the root node */
#define CONSHDLR_PROPFREQ            -1 /**< frequency for propagating domains; zero means only preprocessing propagation */
#define CONSHDLR_EAGERFREQ            1 /**< frequency for using all instead of only the useful constraints in separation,
                                         *   propagation and enforcement, -1 for no eager evaluations, 0 for first only */
#define CONSHDLR_MAXPREROUNDS         0 /**< maximal number of presolving rounds the constraint handler participates in (-1: no limit) */
#define CONSHDLR_DELAYSEPA        FALSE /**< should separation method be delayed, if other separators found cuts? */
#define CONSHDLR_DELAYPROP        FALSE /**< should propagation method be delayed, if other propagators found reductions? */
#define CONSHDLR_DELAYPRESOL      FALSE /**< should presolving method be delayed, if other presolvers found reductions? */
#define CONSHDLR_NEEDSCONS         TRUE /**< should the constraint handler be skipped, if no constraints are available? */




/*
 * Data structures
 */

/** constraint handler data */
struct SCIP_ConshdlrData
{
   /* ???? todo: this should be moved to an appropriat place and has to be assigned in a consistent way (scipex/def.h oder defex.h) */
   mpq_t                 posinfinity;        /**< value considered to be infinity */
   mpq_t                 neginfinity;        /**< value considered to be infinity */
   SCIP_LPIEX*           lpiex;              /**< Exact LP solver interface */
   SCIP_Bool             lpexconstructed;    /**< was the exact LP of some prior node already constructed (constraints)? */
   SCIP_Bool             focuslpexupdated;   /**< was the exact LP of the focus node already updated (boundchanges)?*/
};

/** constraint data for exactlp constraints */
struct SCIP_ConsData
{
   SCIP_OBJSEN           objsense;           /**< objective sense */
   int                   nvars;              /**< number of variables */
   mpq_t*                obj;                /**< objective function values of variables */
   mpq_t*                lb;                 /**< lower bounds of variables */
   mpq_t*                ub;                 /**< upper bounds of variables */
   int                   nconss;             /**< number of constraints */
   mpq_t*                lhs;                /**< left hand sides of constraints */
   mpq_t*                rhs;                /**< right hand sides of constraints */
   int                   nnonz;              /**< number of nonzero elements in the constraint matrix */
   int*                  beg;                /**< start index of each variable in ind- and val-array */
   int*                  ind;                /**< variable indices (var->probindex) of constraint matrix entries */
   mpq_t*                val;                /**< values of nonzero constraint matrix entries */
   SCIP_ROW**            rows;               /**< rows for LP relaxation/approximation (FP data) of exactlp constraint */
};





/*
 * Local methods
 */

/** returns value treated as negative infinite in exactlp constraint handler */
const mpq_t* negInfinity(
   SCIP_CONSHDLRDATA*    conshdlrdata        /**< exactlp constraint handler data */
   )
{  /*lint --e{715} */
   return (const mpq_t*) (&conshdlrdata->neginfinity);
}

/** checks if value is treated as positive infinite in exactlp constraint handler */
const mpq_t* posInfinity(
   SCIP_CONSHDLRDATA*    conshdlrdata        /**< exactlp constraint handler data */
   )
{  /*lint --e{715} */
   return (const mpq_t*) (&conshdlrdata->posinfinity);
}

/** checks if value is treated as negative infinite in exactlp constraint handler */
SCIP_Bool isNegInfinity(
   SCIP_CONSHDLRDATA*    conshdlrdata,       /**< exactlp constraint handler data */
   mpq_t                 val                 /**< value to be compared against infinity */
   )
{  /*lint --e{715} */
   return SCIPlpiexIsNegInfinity(conshdlrdata->lpiex, val);
}

/** checks if value is treated as positive infinite in exactlp constraint handler */
SCIP_Bool isPosInfinity(
   SCIP_CONSHDLRDATA*    conshdlrdata,       /**< exactlp constraint handler data */
   mpq_t                 val                 /**< value to be compared against infinity */
   )
{  /*lint --e{715} */
   return SCIPlpiexIsPosInfinity(conshdlrdata->lpiex, val);
}

/** converts given rational number into an FP number; uses given rounding mode during conversion 
 * (should be used to construct an FP relaxation of a constraint) 
 */
SCIP_Real mpqGetRealRelax(
   SCIP*                 scip,               /**< SCIP data structure */
   mpq_t                 val,                /**< given rational number */
   mp_rnd_t              roundmode           /**< rounding mode to be used for the conversion */
   )
{
   SCIP_Real valrelax;

   mpfr_t valmpfr;
   mpfr_init_set_q(valmpfr, val, roundmode);
   valrelax = (SCIP_Real) mpfr_get_d(valmpfr, roundmode);
   mpfr_clear(valmpfr);

#ifndef NDEBUG
   {
      mpq_t result;
    
      mpq_init(result);
      mpq_set_d(result, valrelax);

      if( roundmode == GMP_RNDU )
         assert(mpq_cmp(result, val) >= 0);
      if( roundmode == GMP_RNDD )
         assert(mpq_cmp(result, val) <= 0); 
#ifdef READER_OUT  /*???????????????*/
      if( roundmode == GMP_RNDU )
      {
         gmp_printf("given <%Qd> <=! conv <%Qd | ", val, result); 
         printf(" %g>\n", mpq_get_d(result)); 
      }
      if( roundmode == GMP_RNDD )
      {
         gmp_printf("given <%Qd> >=! conv <%Qd | ", val, result); 
         printf(" %g>\n", mpq_get_d(result)); 
      }
#endif
   }
#endif

   /* todo: check whether this way to treat infinity is ok (in particular if we want to construct a relaxation) */
   if( SCIPisInfinity(scip, valrelax) )
      valrelax = SCIPinfinity(scip);

   if( SCIPisInfinity(scip, -valrelax) )
      valrelax = -SCIPinfinity(scip);

   return valrelax;
}

/** converts given rational number into an FP number; uses default rounding mode during conversion 
 * (should be used to construct an FP approximation of a constraint) 
 */
SCIP_Real mpqGetRealApprox(
   SCIP*                 scip,               /**< SCIP data structure */
   mpq_t                 val                 /**< given rational number */
   )
{
   SCIP_Real valapprox;

   valapprox = (SCIP_Real) mpq_get_d(val);
 
   /* todo: check whether this way to treat infinity is ok (in particular if we want to construct a relaxation) */
   if( SCIPisInfinity(scip, valapprox) )
      valapprox = SCIPinfinity(scip);

   if( SCIPisInfinity(scip, -valapprox) )
      valapprox = -SCIPinfinity(scip);

   return valapprox;
}

/** creates constaint handler data for exactlp constraint handler */
static
SCIP_RETCODE conshdlrdataCreate(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLRDATA**   conshdlrdata        /**< pointer to store the constraint handler data */
   )
{
   assert(conshdlrdata != NULL);

   SCIP_CALL( SCIPallocMemory(scip, conshdlrdata) );

   QSexactStart(); /* todo: find a good way/place to call this method ?????????????? */

   /* open exact LP Solver interface */
   SCIP_CALL( SCIPlpiexCreate(&(*conshdlrdata)->lpiex, NULL, SCIP_OBJSEN_MINIMIZE) );

   mpq_init((*conshdlrdata)->posinfinity);
   mpq_init((*conshdlrdata)->neginfinity);

   mpq_set((*conshdlrdata)->posinfinity, *SCIPlpiexPosInfinity((*conshdlrdata)->lpiex));
   mpq_set((*conshdlrdata)->neginfinity, *SCIPlpiexNegInfinity((*conshdlrdata)->lpiex));

   (*conshdlrdata)->lpexconstructed = FALSE;
   (*conshdlrdata)->focuslpexupdated = FALSE;

   return SCIP_OKAY;
}

/** frees constraint handler data for linear constraint handler */
static
SCIP_RETCODE conshdlrdataFree(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLRDATA**   conshdlrdata        /**< pointer to the constraint handler data */
   )
{
   assert(conshdlrdata != NULL);
   assert(*conshdlrdata != NULL);

   if( (*conshdlrdata)->lpiex != NULL )
   {
      SCIP_CALL( SCIPlpiexFree(&(*conshdlrdata)->lpiex) );
   }
   assert((*conshdlrdata)->lpiex == NULL);

   mpq_clear((*conshdlrdata)->posinfinity);
   mpq_clear((*conshdlrdata)->neginfinity);

   SCIPfreeMemory(scip, conshdlrdata);

   return SCIP_OKAY;
}

/** gets number of LP rows needed for the LP relaxation of the exactlp constraint */
static
int consdataGetNRows(
   SCIP_CONSDATA*        consdata            /**< constraint data */
   )
{
   assert(consdata != NULL);

   return consdata->nconss;
}

/** creates exactlp constraint data */
static
SCIP_RETCODE consdataCreate(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSDATA**       consdata,           /**< pointer to store constraint data */
   SCIP_OBJSEN           objsense,           /**< objective sense */
   int                   nvars,              /**< number of variables */
   mpq_t*                obj,                /**< objective function values of variables */
   mpq_t*                lb,                 /**< lower bounds of variables */
   mpq_t*                ub,                 /**< upper bounds of variables */
   int                   nconss,              /**< number of constraints */
   mpq_t*                lhs,                /**< left hand sides of constraints */
   mpq_t*                rhs,                /**< right hand sides of constraints */
   int                   nnonz,              /**< number of nonzero elements in the constraint matrix */
   int*                  beg,                /**< start index of each variable in ind- and val-array */
   int*                  ind,                /**< variable indices (var->probindex) of constraint matrix entries */
   mpq_t*                val                 /**< values of nonzero constraint matrix entries */
   )
{
   assert(consdata != NULL);
   assert(nvars > 0 || (nconss == 0 && nconss == 0) );
   assert(nconss > 0 || (nvars >= 0 && nnonz == 0) );
   assert(nnonz > 0 || (nconss >= 0 && nvars >= 0) );

   SCIP_CALL( SCIPallocBlockMemory(scip, consdata) );
   
   /* store variable specific information */ 
   if( nvars > 0 )
   {
      int j;

      /* allocate, initialize, and copy rational arrays */
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(*consdata)->obj, nvars) );
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(*consdata)->lb, nvars) );
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(*consdata)->ub, nvars) );

      for( j = 0; j < nvars; ++j )
      {
         mpq_init((*consdata)->obj[j]);
         mpq_init((*consdata)->lb[j]);
         mpq_init((*consdata)->ub[j]);

         mpq_set((*consdata)->obj[j], obj[j]);
         mpq_set((*consdata)->lb[j], lb[j]);
         mpq_set((*consdata)->ub[j], ub[j]);
      }
   }
   else
   {
      (*consdata)->lb = NULL;
      (*consdata)->ub = NULL;
      (*consdata)->obj = NULL;
   }

   /* store constraint specific information */ 
   if( nconss > 0 )
   {
      int j;

      /* allocate and copy integer array */
      SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &(*consdata)->beg, beg, nconss+1) );

      /* allocate, initialize, and copy rational arrays */
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(*consdata)->lhs, nconss) );
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(*consdata)->rhs, nconss) );

      for( j = 0; j < nconss; ++j )
      {
         mpq_init((*consdata)->lhs[j]);
         mpq_init((*consdata)->rhs[j]);

         mpq_set((*consdata)->lhs[j], lhs[j]);
         mpq_set((*consdata)->rhs[j], rhs[j]);
      }
   }
   else
   {
      (*consdata)->lhs = NULL;
      (*consdata)->rhs = NULL;
      (*consdata)->beg = NULL;
   }

   /* store matrix specific information */ 
   if( nnonz > 0 )
   {
      int j;

      /* allocate and copy integer array */
      SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &(*consdata)->ind, ind, nnonz) );

      /* allocate, initialize, and copy rational array */
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(*consdata)->val, nnonz) );

      for( j = 0; j < nnonz; ++j )
      {
         mpq_init((*consdata)->val[j]);

         mpq_set((*consdata)->val[j], val[j]);
      }
   }
   else
   {
      (*consdata)->ind = NULL;
      (*consdata)->val = NULL;
   }

   (*consdata)->objsense = objsense;
   (*consdata)->nvars = nvars;
   (*consdata)->nconss = nconss;
   (*consdata)->nnonz = nnonz;
   (*consdata)->rows = NULL;

   return SCIP_OKAY;
}

/** releases LP rows of constraint data and frees rows array */
static
SCIP_RETCODE consdataFreeRows(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSDATA*        consdata            /**< constraint data */
   )
{
   assert(consdata != NULL);

   if( consdata->rows != NULL )
   {
      int nrows;
      int r;
      
      nrows = consdataGetNRows(consdata);

      for( r = nrows-1; r >= 0; --r )
      {
         SCIP_CALL( SCIPreleaseRow(scip, &consdata->rows[r]) );
      }
      SCIPfreeBlockMemoryArray(scip, &consdata->rows, nrows);
   }

   return SCIP_OKAY;
}

/** frees exactlp constraint data */
static
SCIP_RETCODE consdataFree(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSDATA**       consdata            /**< pointer to the constraint data */
   )
{
   int j;

   assert(consdata != NULL);
   assert(*consdata != NULL);

   /* release and free the rows */
   SCIP_CALL( consdataFreeRows(scip, *consdata) );

   /* free matrix specific information */ 
   for( j = 0; j < (*consdata)->nnonz; ++j )
   {
      mpq_clear((*consdata)->val[j]);
   }
   SCIPfreeBlockMemoryArrayNull(scip, &(*consdata)->val, (*consdata)->nnonz);
   SCIPfreeBlockMemoryArrayNull(scip, &(*consdata)->ind, (*consdata)->nnonz);

   /* free constraint specific information */ 
   for( j = 0; j < (*consdata)->nconss; ++j )
   {
      mpq_clear((*consdata)->rhs[j]);
      mpq_clear((*consdata)->lhs[j]);
   }
   SCIPfreeBlockMemoryArrayNull(scip, &(*consdata)->rhs, (*consdata)->nconss);
   SCIPfreeBlockMemoryArrayNull(scip, &(*consdata)->lhs, (*consdata)->nconss);
   SCIPfreeBlockMemoryArrayNull(scip, &(*consdata)->beg, (*consdata)->nconss+1);

   /* free variable specific information */ 
   for( j = 0; j < (*consdata)->nvars; ++j )
   {
      mpq_clear((*consdata)->obj[j]);
      mpq_clear((*consdata)->lb[j]);
      mpq_clear((*consdata)->ub[j]);
   }
   SCIPfreeBlockMemoryArrayNull(scip, &(*consdata)->ub, (*consdata)->nvars);
   SCIPfreeBlockMemoryArrayNull(scip, &(*consdata)->lb, (*consdata)->nvars);
   SCIPfreeBlockMemoryArrayNull(scip, &(*consdata)->obj, (*consdata)->nvars);

   /* free consdata */ 
   SCIPfreeBlockMemory(scip, consdata);
 
   return SCIP_OKAY;
}

/** prints exactlp constraint in CIP format to file stream */
static
void consdataPrint(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLRDATA*    conshdlrdata,       /**< exactlp constraint handler data */
   SCIP_CONSDATA*        consdata,           /**< linear constraint data */
   FILE*                 file                /**< output file (or NULL for standard output) */
   )
{
   SCIP_VAR** vars;
   char s[SCIP_MAXSTRLEN];
   int i;

   assert(conshdlrdata != NULL);
   assert(consdata != NULL);

   /* get problem variables */
   vars = SCIPgetOrigVars(scip);
   assert(SCIPgetNOrigVars(scip) == consdata->nvars);

   SCIPinfoMessage(scip, file, "\n");

   /* print objective sense */
   SCIPinfoMessage(scip, file, "  OBJECTIVE EXACT\n");
   SCIPinfoMessage(scip, file, "    Sense            : %s\n", consdata->objsense == SCIP_OBJSENSE_MINIMIZE ? "minimize" : "maximize");
   
   /* print variable specific information */
   if( consdata->nvars > 0 )
   {
      SCIPinfoMessage(scip, file, "  VARIABLES EXACT\n");
      for( i = 0; i < consdata->nvars; ++i )
      {
         assert(SCIPvarGetProbindex(vars[i]) == i);

         /* name */
         SCIPmessageFPrintInfo(file, "    <%s>:", SCIPvarGetName(vars[i]));
         
         /* objective value */
         gmp_snprintf(s, SCIP_MAXSTRLEN, " obj=%Qd", consdata->obj[i]);
         SCIPinfoMessage(scip, file, s);

         /* bounds (global bounds for transformed variables, original bounds for original variables) */
         SCIPmessageFPrintInfo(file, ", bounds=");

         if( isPosInfinity(conshdlrdata, consdata->lb[i]) )
            SCIPmessageFPrintInfo(file, "[+inf,");
         else if( isNegInfinity(conshdlrdata, consdata->lb[i]) )
            SCIPmessageFPrintInfo(file, "[-inf,");
         else
         {
            gmp_snprintf(s, SCIP_MAXSTRLEN, "[%Qd,", consdata->lb[i]);
            SCIPinfoMessage(scip, file, s);
         }

         if( isPosInfinity(conshdlrdata, consdata->ub[i]) )
            SCIPmessageFPrintInfo(file, "+inf]");
         else if( isNegInfinity(conshdlrdata, consdata->ub[i]) )
            SCIPmessageFPrintInfo(file, "-inf]");
         else
         {
            gmp_snprintf(s, SCIP_MAXSTRLEN, "%Qd]", consdata->ub[i]);
            SCIPinfoMessage(scip, file, s);
         }

         SCIPmessageFPrintInfo(file, "\n");
      }
   }

   /* print constraint and matrix specific information */
   if( consdata->nconss > 0 )
   {
      SCIPinfoMessage(scip, file, "  CONSTRAINTS EXACT\n");

      for( i = 0; i < consdata->nconss; ++i )
      {
         SCIPinfoMessage(scip, file, "    ");

         /* print left hand side for ranged rows */
         if( !isNegInfinity(conshdlrdata, consdata->lhs[i])
            && !isPosInfinity(conshdlrdata, consdata->rhs[i])
            && mpq_equal(consdata->lhs[i], consdata->rhs[i]) == 0 )
         {
            gmp_snprintf(s, SCIP_MAXSTRLEN, "%Qd <= ", consdata->lhs[i]);
            SCIPinfoMessage(scip, file, s);
         }
         
         /* print coefficients and variables */
         if( consdata->beg[i] == consdata->beg[i+1] )
            SCIPinfoMessage(scip, file, "0 ");
         else
         {
            int v;
            
            for( v = consdata->beg[i]; v < consdata->beg[i+1]; v++ )
            {
               assert(consdata->ind[v] >= 0 && consdata->ind[v] < consdata->nvars);
               assert(SCIPvarGetProbindex(vars[consdata->ind[v]]) == consdata->ind[v]);
 
               gmp_snprintf(s, SCIP_MAXSTRLEN, "%+Qd<%s> ", consdata->val[v], SCIPvarGetName(vars[consdata->ind[v]]));
               SCIPinfoMessage(scip, file, s);
            }
         }

         /* print right hand side */
         if( mpq_equal(consdata->lhs[i], consdata->rhs[i]) != 0 )
         {
            gmp_snprintf(s, SCIP_MAXSTRLEN, "== %Qd\n", consdata->rhs[i]);
            SCIPinfoMessage(scip, file, s);
         }
         else if( !isPosInfinity(conshdlrdata, consdata->rhs[i]) )
         {
            gmp_snprintf(s, SCIP_MAXSTRLEN, "<= %Qd\n", consdata->rhs[i]);
            SCIPinfoMessage(scip, file, s);
         }
         else if( !isNegInfinity(conshdlrdata, consdata->lhs[i]) )
         {
            gmp_snprintf(s, SCIP_MAXSTRLEN, ">= %Qd\n", consdata->lhs[i]);
            SCIPinfoMessage(scip, file, s);
         }
         else
            SCIPinfoMessage(scip, file, " [free]\n");
      }
   }
}

/** creates LP rows corresponding to exactlp constraint */
static 
SCIP_RETCODE createRelaxation(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint to check */
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONSDATA* consdata;
   SCIP_VAR** vars;
   SCIP_VAR** rowvars; 
   SCIP_Real* rowvals; 
   SCIP_Real rowlhs;
   SCIP_Real rowrhs;
   int nrowvars;
   
   char rowname[SCIP_MAXSTRLEN];
   int nrows;
   int c;
   int i;

   conshdlrdata = SCIPconshdlrGetData(SCIPconsGetHdlr(cons));
   assert(conshdlrdata != NULL);
   
   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);
   assert(consdata->rows == NULL);

   /* get problem variables */
   vars = SCIPgetVars(scip);
   nrows = consdataGetNRows(consdata);
   assert(SCIPgetNVars(scip) == consdata->nvars);
   assert(nrows == consdata->nconss);

   /* allocate memory for all rows */
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &consdata->rows, nrows) );

   /* allocate temporary memory */
   SCIP_CALL( SCIPallocBufferArray(scip, &rowvars, SCIPgetNVars(scip)) );
   SCIP_CALL( SCIPallocBufferArray(scip, &rowvals, SCIPgetNVars(scip)) );
   
#ifdef CONSEXLP_OUT /* only for debugging ???????????????? */
   SCIP_CALL( SCIPprintCons(scip, cons, NULL) ); 
#endif

   if( SCIPuseFPRelaxation(scip) ) 
   {
     /* for each row of the exactlp constraint, create a row with FP data that defines a relaxation */
      for( c = 0; c < consdata->nconss; ++c )
      {
         int v;

         /* calculate lhs and rhs of row */
         rowlhs = mpqGetRealRelax(scip, consdata->lhs[c], GMP_RNDD);
         rowrhs = mpqGetRealRelax(scip, consdata->rhs[c], GMP_RNDU);

         if( !SCIPisInfinity(scip, -rowlhs) && !SCIPisInfinity(scip, rowrhs) )
         {
            /* todo: create two rows and modify consdata structures in order to be able to transform a basis 
             * to the original problem later ??????? 
             */
            SCIPerrorMessage("consinitlp: for ranged rows, creating a FP relaxation is not supported yet\n");
            SCIPABORT(); /*lint --e{527}*/
            goto TERMINATE;
         }            
        
         v = 0;
         nrowvars = consdata->beg[c+1] - consdata->beg[c];
         assert(nrowvars >= 0 && nrowvars <= SCIPgetNVars(scip));

         /* calculate coefficients of all variable in the row */
         for( i = consdata->beg[c]; i < consdata->beg[c+1]; ++i )
         {
            int probidx;

            probidx = consdata->ind[i];
            
            assert(SCIPvarGetProbindex(vars[probidx]) == probidx);
            assert(consdata->objsense * mpqGetRealApprox(scip, consdata->obj[probidx]) == SCIPvarGetObj(vars[probidx]));
            assert(mpqGetRealRelax(scip, consdata->lb[probidx], GMP_RNDD) == SCIPvarGetLbGlobal(vars[probidx]));
            assert(mpqGetRealRelax(scip, consdata->ub[probidx], GMP_RNDU) == SCIPvarGetUbGlobal(vars[probidx]));
            assert(mpqGetRealRelax(scip, consdata->lb[probidx], GMP_RNDD) == SCIPvarGetLbLocal(vars[probidx]));
            assert(mpqGetRealRelax(scip, consdata->ub[probidx], GMP_RNDU) == SCIPvarGetUbLocal(vars[probidx]));
            assert(mpq_cmp(consdata->lb[probidx], consdata->ub[probidx]) <= 0);

            /* x_j >= 0 holds always ==> 
             *   cons with rhs: underestimate(a_j) * x_j <= a_j * x_j <= rhs
             *   cons with lhs:  overestimate(a_j) * x_j >= a_j * x_j >= lhs 
             */
            if( mpq_sgn(consdata->lb[probidx]) >= 0 )
            {
               if( !SCIPisInfinity(scip, rowrhs) )
               {
                  assert(SCIPisInfinity(scip, -rowlhs));
                  rowvals[v] = mpqGetRealRelax(scip, consdata->val[i], GMP_RNDD);
               }
               else
               {
                  assert(!SCIPisInfinity(scip, -rowlhs) && SCIPisInfinity(scip, rowrhs));
                  rowvals[v] = mpqGetRealRelax(scip, consdata->val[i], GMP_RNDU);
               }
            }
            /* x_j <= 0 always holds ==>
             *   cons with rhs:  overestimate(a_j) * x_j <= a_j * x_j <= rhs
             *   cons with lhs: underestimate(a_j) * x_j >= a_j * x_j >= lhs
             */
            else if( mpq_sgn(consdata->ub[probidx]) <= 0 )
            {
               if( !SCIPisInfinity(scip, rowrhs) )
               {
                  assert(SCIPisInfinity(scip, -rowlhs));
                  rowvals[v] = mpqGetRealRelax(scip, consdata->val[i], GMP_RNDU);
               }
               else
               {
                  assert(!SCIPisInfinity(scip, -rowlhs) && SCIPisInfinity(scip, rowrhs));
                  rowvals[v] = mpqGetRealRelax(scip, consdata->val[i], GMP_RNDD);
               }
            }
            /* x_j <= 0 and x_j >= 0 may hold ==> 
             *   split x_j into negative and positive part 
             */
            else
            {
               /* todo: split variable into positive and negative part and modify consdata structures in order to be able 
                * to transform a basis to the original problem later ??????? 
                */
               SCIPerrorMessage("consinitlp: for variable that are neither nonnegative nor nonpositive, creating a FP relaxation is not supported yet\n");
               SCIPABORT(); /*lint --e{527}*/
               goto TERMINATE;
            }

            rowvars[v] = vars[probidx];
            v++;
         }
         
         /* create the row */
         (void) SCIPsnprintf(rowname, SCIP_MAXSTRLEN, "%s_relax_%d", SCIPconsGetName(cons), c);
         SCIP_CALL( SCIPcreateEmptyRow(scip, &consdata->rows[c], rowname, rowlhs, rowrhs, SCIPconsIsLocal(cons), 
               SCIPconsIsModifiable(cons), SCIPconsIsRemovable(cons)) );
         SCIP_CALL( SCIPaddVarsToRow(scip, consdata->rows[c], nrowvars, rowvars, rowvals) );

#ifdef CONSEXLP_OUT /* only for debugging ???????????????? */
         SCIPdebug(SCIPprintRow(scip, consdata->rows[c], NULL)); 
#endif
      }
   }
   else
   {
      /* for each row of the exactlp constraint, create a row with FP data that defines an approximation */
      for( c = 0; c < consdata->nconss; ++c )
      {
         int v;

         /* calculate lhs and rhs of row */
         rowlhs = mpqGetRealApprox(scip, consdata->lhs[c]);
         rowrhs = mpqGetRealApprox(scip, consdata->rhs[c]);

         v = 0;
         nrowvars = consdata->beg[c+1] - consdata->beg[c];
         assert(nrowvars >= 0 && nrowvars <= SCIPgetNVars(scip));

         /* add all variables to the row */
         for( i = consdata->beg[c]; i < consdata->beg[c+1]; ++i )
         {
            int probidx;
            probidx = consdata->ind[i];

            assert(SCIPvarGetProbindex(vars[probidx]) == probidx);
            assert(consdata->objsense * mpqGetRealApprox(scip, consdata->obj[probidx]) == SCIPvarGetObj(vars[probidx]));
            assert(mpqGetRealRelax(scip, consdata->lb[probidx], GMP_RNDD) == SCIPvarGetLbGlobal(vars[probidx]));
            assert(mpqGetRealRelax(scip, consdata->ub[probidx], GMP_RNDU) == SCIPvarGetUbGlobal(vars[probidx]));
            assert(mpqGetRealRelax(scip, consdata->lb[probidx], GMP_RNDD) == SCIPvarGetLbLocal(vars[probidx]));
            assert(mpqGetRealRelax(scip, consdata->ub[probidx], GMP_RNDU) == SCIPvarGetUbLocal(vars[probidx]));
            
            rowvals[v] = mpqGetRealApprox(scip, consdata->val[i]);
            rowvars[v] = vars[probidx];
            v++;
         }

         /* create the row */
         (void) SCIPsnprintf(rowname, SCIP_MAXSTRLEN, "%s_approx_%d", SCIPconsGetName(cons), c);
         SCIP_CALL( SCIPcreateEmptyRow(scip, &consdata->rows[c], rowname, rowlhs, rowrhs, SCIPconsIsLocal(cons), 
               SCIPconsIsModifiable(cons), SCIPconsIsRemovable(cons)) );
         SCIP_CALL( SCIPaddVarsToRow(scip, consdata->rows[c], nrowvars, rowvars, rowvals) );
         
#ifdef CONSEXLP_OUT /* only for debugging ???????????????? */
         SCIPdebug(SCIPprintRow(scip, consdata->rows[c], NULL)); 
#endif
      }
   }

 TERMINATE:
   /* free temporary memory */
   SCIPfreeBufferArray(scip, &rowvals);
   SCIPfreeBufferArray(scip, &rowvars);

   return SCIP_OKAY;
}  

/** adds linear relaxation of exactlp constraint to the LP */
static 
SCIP_RETCODE addRelaxation(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint to check */
   )
{
   SCIP_CONSDATA* consdata;
   int nrows;
   int r;

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   if( consdata->rows == NULL )
   {
      SCIP_CALL( createRelaxation(scip, cons) );
   }
   assert( consdata->rows != NULL );

   nrows = consdataGetNRows(consdata);

   for( r = 0; r < nrows; ++r )
   {
      if( !SCIProwIsInLP(consdata->rows[r]) )
      {
         /* todo: check whether it is ok/better not to force the cut to enter the LP ???????????? */
         SCIP_CALL( SCIPaddCut(scip, NULL, consdata->rows[r], TRUE) ); 
      }
   }

   return SCIP_OKAY;
}

/** constructs the exact LP of the current node, but does not load the LP state and warmstart information  */
static
SCIP_RETCODE constructCurrentLPEX(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLRDATA*    conshdlrdata,       /**< exactlp constraint handler data */
   SCIP_CONSDATA*        consdata            /**< exactlp constraint data */
   )
{
   assert(conshdlrdata != NULL);
   assert(consdata != NULL);

   SCIPdebugMessage("constructing initial exact LP\n");

   if( !conshdlrdata->lpexconstructed )
   {
      SCIP_VAR** vars;
      int i;
      char** colnames; /* todo: is this implemented in a correct way ????????????*/
       
      /* allocate and initialize temporary memory */
      SCIP_CALL( SCIPallocBufferArray(scip, &colnames, consdata->nvars) );

      /* get names of problem variables */
      vars = SCIPgetVars(scip);
      assert(SCIPgetNVars(scip) == consdata->nvars);
      for( i = 0; i < consdata->nvars; ++i )
      {
         /* allocate and initialize temporary memory */
         colnames[i] = (char*) (SCIPvarGetName(vars[i]));

         assert(SCIPvarGetProbindex(vars[i]) == i);
         assert(consdata->objsense * mpqGetRealApprox(scip, consdata->obj[i]) == SCIPvarGetObj(vars[i]));
         assert(mpqGetRealRelax(scip, consdata->lb[i], GMP_RNDD) <= SCIPvarGetLbGlobal(vars[i]));
         assert(mpqGetRealRelax(scip, consdata->ub[i], GMP_RNDU) >= SCIPvarGetUbGlobal(vars[i]));
         assert(mpqGetRealRelax(scip, consdata->lb[i], GMP_RNDD) <= SCIPvarGetLbLocal(vars[i]));
         assert(mpqGetRealRelax(scip, consdata->ub[i], GMP_RNDU) >= SCIPvarGetUbLocal(vars[i]));
      }

      /* add all columns to the exact LP */
      /* todo: check whether I implement SCIPlpiexAddCols() correctly (handling of case: beg=ind=val=NULL) ????????????? */
      SCIP_CALL( SCIPlpiexAddCols(conshdlrdata->lpiex, consdata->nvars, consdata->obj, consdata->lb, consdata->ub, colnames,
            0, NULL, NULL, NULL) );

      /* add all constraints to the exact LP */
      SCIP_CALL( SCIPlpiexAddRows(conshdlrdata->lpiex, consdata->nconss, consdata->lhs, consdata->rhs, NULL, 
            consdata->nnonz, consdata->beg, consdata->ind, consdata->val) );

      // SCIP_CALL( SCIPlpiexWriteLP(conshdlrdata->lpiex, "testset/debug2.lp") ); /* ????????????*/
      // SCIPABORT(); /*lint --e{527}*/ /* ????????????*/

      conshdlrdata->lpexconstructed = TRUE;

      /* free temporary memory */
      SCIPfreeBufferArray(scip, &colnames);
   }

   /* todo: CONTINUE here (later): ????????????
    * - implement eventhandler for boundchanges and update bounds of variables 
    *   (decide whether to update consdata->lb/ub arrays in addition)
    */
   if( !conshdlrdata->focuslpexupdated )
   {
      SCIPerrorMessage("constructCurrentLPEX: updateboundchanges part of method of exactlp constraint handler not implemented yet\n");
      SCIPABORT(); /*lint --e{527}*/
   }

   return SCIP_OKAY;
}

/** solves the exact LP with the given algorithm and evaluates return status */
static
SCIP_RETCODE solveLPEX(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLRDATA*    conshdlrdata,       /**< exactlp constraint handler data */
   SCIP_LPALGO           lpalgo,             /**< LP algorithm that should be applied */
   SCIP_Bool*            lperror             /**< pointer to store whether an unresolved LP error occured */
   )
{
   SCIP_RETCODE retcode;
   int iterations;
   int ncols;
   int nrows;

   assert(scip != NULL);
   assert(lperror != NULL);

   
   SCIPlpiexGetNCols(conshdlrdata->lpiex, &ncols);
   SCIPlpiexGetNRows(conshdlrdata->lpiex, &nrows);
   
   *lperror = FALSE;

   switch( lpalgo )
   {
   case SCIP_LPALGO_PRIMALSIMPLEX:
      SCIPdebugMessage("solving current primal exact LP (%d cols, %d rows)\n", ncols, nrows);

      /* call primal simplex */
      retcode = SCIPlpiexSolvePrimal(conshdlrdata->lpiex);
      if( retcode == SCIP_LPERROR )
      {
         *lperror = TRUE;
         SCIPdebugMessage("(node %"SCIP_LONGINT_FORMAT") primal simplex solving error in current exact LP\n", 
            SCIPgetNNodes(scip));
      }
      else
      {
         SCIP_CALL( retcode );
      }
 
      SCIP_CALL( SCIPlpiexGetIterations(conshdlrdata->lpiex, &iterations) );
      SCIPdebugMessage("solved primal exact LP in %d iterations\n", iterations);
      break;

   case SCIP_LPALGO_DUALSIMPLEX:
      SCIPdebugMessage("solving current dual exact LP (%d cols, %d rows)\n", ncols, nrows);
      
      /* call dual simplex */
      retcode = SCIPlpiexSolveDual(conshdlrdata->lpiex);
      if( retcode == SCIP_LPERROR )
      {
         *lperror = TRUE;
         SCIPdebugMessage("(node %"SCIP_LONGINT_FORMAT") dual simplex solving error in current exact LP\n", 
            SCIPgetNNodes(scip));
      }
      else
      {
         SCIP_CALL( retcode );
      }

      SCIP_CALL( SCIPlpiexGetIterations(conshdlrdata->lpiex, &iterations) );
      SCIPdebugMessage("solved dual exact LP in %d iterations\n", iterations);
      break;

   default:
      SCIPerrorMessage("invalid exact LP algorithm\n");
      return SCIP_INVALIDDATA;
   }

   return SCIP_OKAY;
}

/*
 * Callback methods of constraint handler
 */

/** destructor of constraint handler to free constraint handler data (called when SCIP is exiting) */
static
SCIP_DECL_CONSFREE(consFreeExactlp)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;

   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);

   /* free constraint handler data */
   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   SCIP_CALL( conshdlrdataFree(scip, &conshdlrdata) );

   SCIPconshdlrSetData(conshdlr, NULL);

   return SCIP_OKAY;
}


/** initialization method of constraint handler (called after problem was transformed) */
#define consInitExactlp NULL


/** deinitialization method of constraint handler (called before transformed problem is freed) */
#define consExitExactlp NULL


/** presolving initialization method of constraint handler (called when presolving is about to begin) */
#if 0
static
SCIP_DECL_CONSINITPRE(consInitpreExactlp)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of exactlp constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consInitpreExactlp NULL
#endif


/** presolving deinitialization method of constraint handler (called after presolving has been finished) */
#if 0
static
SCIP_DECL_CONSEXITPRE(consExitpreExactlp)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of exactlp constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consExitpreExactlp NULL
#endif


/** solving process initialization method of constraint handler (called when branch and bound process is about to begin) */
#if 0
static
SCIP_DECL_CONSINITSOL(consInitsolExactlp)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of exactlp constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consInitsolExactlp NULL
#endif


/** solving process deinitialization method of constraint handler (called before branch and bound process data is freed) */
#if 0
static
SCIP_DECL_CONSEXITSOL(consExitsolExactlp)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of exactlp constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consExitsolExactlp NULL
#endif


/** frees specific constraint data */
static
SCIP_DECL_CONSDELETE(consDeleteExactlp)
{  /*lint --e{715}*/
   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);
 
   /* free exactlp constraint */
   SCIP_CALL( consdataFree(scip, consdata) );
   
   return SCIP_OKAY;
}


/** transforms constraint data into data belonging to the transformed problem */ 
static
SCIP_DECL_CONSTRANS(consTransExactlp)
{  /*lint --e{715}*/
   SCIP_CONSDATA* sourcedata;
   SCIP_CONSDATA* targetdata;

   /*debugMessage("Trans method of exactlp constraints\n");*/

   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);
   assert(SCIPgetStage(scip) == SCIP_STAGE_TRANSFORMING);
   assert(sourcecons != NULL);
   assert(targetcons != NULL);

   sourcedata = SCIPconsGetData(sourcecons);
   assert(sourcedata != NULL);

   /* create exactlp constraint data for target constraint */
   SCIP_CALL( consdataCreate(scip, &targetdata, sourcedata->objsense, sourcedata->nvars, sourcedata->obj, sourcedata->lb, 
         sourcedata->ub, sourcedata->nconss, sourcedata->lhs, sourcedata->rhs, sourcedata->nnonz, sourcedata->beg, 
         sourcedata->ind, sourcedata->val) );
 
   /* create target constraint */
   SCIP_CALL( SCIPcreateCons(scip, targetcons, SCIPconsGetName(sourcecons), conshdlr, targetdata,
         SCIPconsIsInitial(sourcecons), SCIPconsIsSeparated(sourcecons), SCIPconsIsEnforced(sourcecons),
         SCIPconsIsChecked(sourcecons), SCIPconsIsPropagated(sourcecons),
         SCIPconsIsLocal(sourcecons), SCIPconsIsModifiable(sourcecons),
         SCIPconsIsDynamic(sourcecons), SCIPconsIsRemovable(sourcecons), SCIPconsIsStickingAtNode(sourcecons)) );

   return SCIP_OKAY;
}


/** LP initialization method of constraint handler */
static
SCIP_DECL_CONSINITLP(consInitlpExactlp)
{  /*lint --e{715}*/
   int i;

   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);

   for( i = 0; i < nconss; ++i )
   {
      assert(SCIPconsIsInitial(conss[i]));
      SCIP_CALL( addRelaxation(scip, conss[i]) );
   }

   return SCIP_OKAY;
}

/** separation method of constraint handler for LP solutions */
/* todo: implement this (because of sepafreq = 1, there will be an error ??????? */
static
SCIP_DECL_CONSSEPALP(consSepalpExactlp)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONSDATA* consdata;

   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);
   assert(result != NULL);
   assert(nconss <= 1);

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   consdata = SCIPconsGetData(conss[0]);
   assert(consdata != NULL);

   *result = SCIP_DIDNOTRUN; /* todo: find out which result should be returend ?????????*/

   /* in case the FP problem is a relaxation of the original problem, we have already calculated a proved lower bound
    * via postprocessing the LP solution of the FP problem 
    */
   if( SCIPuseFPRelaxation(scip) ) 
      return SCIP_OKAY;

   /* todo: CONTINUE here ?????????
    * - fill lpiex
    * - calculate safe lower bound and hand it to SCIP
    */
   //   SCIP_CALL( constructCurrentLPEX(scip, conshdlrdata, consdata) );
   
   return SCIP_OKAY;
}

/** separation method of constraint handler for arbitrary primal solutions */
#if 0
static
SCIP_DECL_CONSSEPASOL(consSepasolExactlp)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of exactlp constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consSepasolExactlp NULL
#endif


/** constraint enforcing method of constraint handler for LP solutions */
static
SCIP_DECL_CONSENFOLP(consEnfolpExactlp)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONSDATA* consdata;
   SCIP_Bool lperror; 
   int ncols;
   char algo;

   assert(SCIPhasCurrentNodeLP(scip));
   assert(SCIPgetLPSolstat(scip) == SCIP_LPSOLSTAT_OPTIMAL);

   SCIPdebugMessage("enforcing exactlp constraint <%s>\n", SCIPconsGetName(conss[0]));

   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);
   assert(result != NULL);
   assert(nconss <= 1);

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   consdata = SCIPconsGetData(conss[0]);
   assert(consdata != NULL);

   *result = SCIP_FEASIBLE;

   conshdlrdata->focuslpexupdated = TRUE; /* todo: just a workaround in order to solve the first lpex ??????? */
   constructCurrentLPEX(scip, conshdlrdata, consdata);

   SCIPlpiexGetNCols(conshdlrdata->lpiex, &ncols);
   assert(SCIPgetNVars(scip) == ncols);
    
   /* todo: implement the enforcement method ????????????
    * - spaeter: warmstartinfo/basis an LPEX uebergeben und 
    *            dual/primal simplex entsprechend das status der uebergebenen basis entscheiden
    */

   /* solve the exact LP relaxation */
   algo = 'd'; /* todo: choose wrt given warmstart basis */
   switch( algo )
   {
   case 'd':
      SCIP_CALL( solveLPEX(scip, conshdlrdata, SCIP_LPALGO_DUALSIMPLEX, &lperror) );
      break;

   case 'p':
      SCIP_CALL( solveLPEX(scip, conshdlrdata, SCIP_LPALGO_PRIMALSIMPLEX, &lperror) );
      break;

   default:
      SCIPerrorMessage("invalid parameter setting <%c> for exact LP algorithm\n", algo);
      return SCIP_PARAMETERWRONGVAL;
   }

   /* evaluate solution status */
   if( lperror )
   {
      SCIPerrorMessage("exact LP solver returns error: case not handled yet\n");
      SCIPABORT();
   }

   if ( SCIPlpiexIsOptimal(conshdlrdata->lpiex) )
   {
      SCIPdebugMessage("exact LP solved to optimality\n"); 

      SCIP_Bool primalfeasible;
      SCIP_Bool dualfeasible;
      mpq_t* primsol;
      mpq_t lpobjval;
      int j;

      /* check for primal and dual feasibility */
      SCIP_CALL( SCIPlpiexGetSolFeasibility(conshdlrdata->lpiex, &primalfeasible, &dualfeasible) );
      assert(primalfeasible);
      assert(dualfeasible);
      
      /* allocate and initialize temporary arrays */
      SCIP_CALL( SCIPallocBufferArray(scip, &primsol, ncols) );
      for( j = 0; j < ncols; ++j )
      {
         mpq_init(primsol[j]);
      }

      /* check whether primal solution satisfies all integrality restrictions */
      SCIP_CALL( SCIPlpiexGetSol(conshdlrdata->lpiex, &lpobjval, primsol, NULL, NULL, NULL) );
      
      /* todo: CONTINUE here: implement the enforcement method ????????????
       * - auf ganzzahligkeit testen (step 3b)
       */
      //      checkIntegrality(scip, primsol);

      /* free temporary arrays */
      for( j = 0; j < ncols; ++j )
      {
         mpq_clear(primsol[j]);
      }
      SCIPfreeBufferArray(scip, &primsol);
   }
   else if( SCIPlpiexIsObjlimExc(conshdlrdata->lpiex) )
   {
      SCIPerrorMessage("exact LP exceeds objlimit: case not handled yet\n");
      SCIPABORT();
   }
   else if( SCIPlpiexIsPrimalInfeasible(conshdlrdata->lpiex) )
   {
      /* todo: CONTINUE here: implement the enforcement method ????????????
       * - knoten abschneiden (step 3c)
       */
      SCIPerrorMessage("exact LP primal infeasible: case not handled yet\n");
      SCIPABORT();
   }
   else if( SCIPlpiexExistsPrimalRay(conshdlrdata->lpiex) ) /* todo: check why in lp.c SCIPlpiIsPrimalUnbounded() is not used ????*/
   {
      SCIPerrorMessage("exact LP has primal ray: case not handled yet\n");
      SCIPABORT();
   }
   else if( SCIPlpiexIsIterlimExc(conshdlrdata->lpiex) )
   {
      SCIPerrorMessage("exact LP exceeds iteration limit: case not handled yet\n");
      SCIPABORT();
   }
   else if( SCIPlpiexIsTimelimExc(conshdlrdata->lpiex) )
   {
      SCIPerrorMessage("exact LP exceeds time limit: case not handled yet\n");
      SCIPABORT();
   }
   else
   {
      SCIPerrorMessage("(node %"SCIP_LONGINT_FORMAT") error or unknown return status in current exact LP (internal status: %d)\n", 
         SCIPgetNNodes(scip), SCIPlpiexGetInternalStatus(conshdlrdata->lpiex));
      return SCIP_LPERROR;
   }

   SCIPABORT();
   *result = SCIP_INFEASIBLE; /* only because methode is not completely implemented yet ??????????? */

   return SCIP_OKAY;
}


/** constraint enforcing method of constraint handler for pseudo solutions */
static
SCIP_DECL_CONSENFOPS(consEnfopsExactlp)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of exactlp constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}


/** feasibility check method of constraint handler for integral solutions */
static
SCIP_DECL_CONSCHECK(consCheckExactlp)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of exactlp constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}


/** domain propagation method of constraint handler */
#if 0
static
SCIP_DECL_CONSPROP(consPropExactlp)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of exactlp constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consPropExactlp NULL
#endif


/** presolving method of constraint handler */
#if 0
static
SCIP_DECL_CONSPRESOL(consPresolExactlp)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of exactlp constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consPresolExactlp NULL
#endif


/** propagation conflict resolving method of constraint handler */
#if 0
static
SCIP_DECL_CONSRESPROP(consRespropExactlp)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of exactlp constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consRespropExactlp NULL
#endif


/** variable rounding lock method of constraint handler */
static
SCIP_DECL_CONSLOCK(consLockExactlp)
{  /*lint --e{715}*/
   SCIP_VAR** vars;
   SCIP_CONSDATA* consdata;
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_Bool* roundinguplocked;
   SCIP_Bool* roundingdownlocked;
   SCIP_Bool haslhs;
   SCIP_Bool hasrhs;
   int j;
   int c;

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   /* get problem variables; in transforming stage we can only access original variables, 
    * however, SCIPaddVarLocks() will use the transformed variable if it exists 
    */
   vars = SCIPgetOrigVars(scip); 
   assert(SCIPgetNOrigVars(scip) == consdata->nvars);

   /* allocate and initialize temporary memory */
   SCIP_CALL( SCIPallocBufferArray(scip, &roundingdownlocked, consdata->nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &roundinguplocked, consdata->nvars) );
   BMSclearMemoryArray(roundingdownlocked, consdata->nvars);
   BMSclearMemoryArray(roundinguplocked, consdata->nvars);

   /* for every variable, check whether rounding up or down could lead to the violation of some constraint */
   j = 0; 
   for( c = 0; c < consdata->nconss; ++c )
   {
      haslhs = !isNegInfinity(conshdlrdata, consdata->lhs[c]);
      hasrhs = !isPosInfinity(conshdlrdata, consdata->rhs[c]);

      /* check all variables of the current constraint */
      for( ; j < consdata->beg[c+1]; ++j )
      {
         assert(consdata->ind[j] >= 0 && consdata->ind[j] < consdata->nvars);

         /* val > 0 */
         if( mpq_sgn(consdata->val[j]) > 0 )
         {
            if( haslhs )
               roundingdownlocked[consdata->ind[j]] = TRUE;
            
            if( hasrhs )
               roundinguplocked[consdata->ind[j]] = TRUE;
         }

         /* val < 0 */
         if( mpq_sgn(consdata->val[j]) < 0 )
         {
            if( haslhs )
               roundinguplocked[consdata->ind[j]] = TRUE;
            
            if( hasrhs )
               roundingdownlocked[consdata->ind[j]] = TRUE;
         }
      } 
   }
   assert(j == consdata->nnonz);

   /* set rounding locks for all variables */
   for( j = 0; j < consdata->nvars; ++j )
   {
      int probindex;

      probindex = SCIPvarGetProbindex(vars[j]);

      if( roundingdownlocked[probindex] && roundinguplocked[probindex] )
      {
         SCIP_CALL( SCIPaddVarLocks(scip, vars[j], nlockspos + nlocksneg, nlockspos + nlocksneg) );
      }
      else if( roundingdownlocked[probindex] )
      { 
         assert(!roundinguplocked[probindex]);

         SCIP_CALL( SCIPaddVarLocks(scip, vars[j], nlockspos, nlocksneg) );
      } 
      else if( roundinguplocked[probindex] )
      { 
         assert(!roundingdownlocked[probindex]);

         SCIP_CALL( SCIPaddVarLocks(scip, vars[j], nlocksneg, nlockspos) );
      } 
   }
   
   /* free temporary memory */
   SCIPfreeBufferArray(scip, &roundinguplocked);
   SCIPfreeBufferArray(scip, &roundingdownlocked);

   return SCIP_OKAY;
}


/** constraint activation notification method of constraint handler */
#if 0
static
SCIP_DECL_CONSACTIVE(consActiveExactlp)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of exactlp constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consActiveExactlp NULL
#endif


/** constraint deactivation notification method of constraint handler */
#if 0
static
SCIP_DECL_CONSDEACTIVE(consDeactiveExactlp)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of exactlp constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consDeactiveExactlp NULL
#endif


/** constraint enabling notification method of constraint handler */
#if 0
static
SCIP_DECL_CONSENABLE(consEnableExactlp)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of exactlp constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consEnableExactlp NULL
#endif


/** constraint disabling notification method of constraint handler */
#if 0
static
SCIP_DECL_CONSDISABLE(consDisableExactlp)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of exactlp constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consDisableExactlp NULL
#endif


/** constraint display method of constraint handler */
static
SCIP_DECL_CONSPRINT(consPrintExactlp)
{  /*lint --e{715}*/
   assert( scip != NULL );
   assert( conshdlr != NULL );
   assert( cons != NULL );
   
   consdataPrint(scip, SCIPconshdlrGetData(conshdlr), SCIPconsGetData(cons), file);
    
   return SCIP_OKAY;
}




/*
 * Linear constraint upgrading
 */

#ifdef LINCONSUPGD_PRIORITY
/** tries to upgrade a linear constraint into a exactlp constraint */
static
SCIP_DECL_LINCONSUPGD(linconsUpgdExactlp)
{  /*lint --e{715}*/
   SCIP_Bool upgrade;

   assert(upgdcons != NULL);
   
   /* check, if linear constraint can be upgraded to exactlp constraint */
   upgrade = FALSE;
   /* TODO: put the constraint's properties here, in terms of the statistics given by nposbin, nnegbin, ... */

   if( upgrade )
   {
      SCIPdebugMessage("upgrading constraint <%s> to exactlp constraint\n", SCIPconsGetName(cons));
      
      /* create the bin Exactlp constraint (an automatically upgraded constraint is always unmodifiable) */
      assert(!SCIPconsIsModifiable(cons));
      SCIP_CALL( SCIPcreateConsExactlp(scip, upgdcons, SCIPconsGetName(cons), nvars, vars, vals, lhs, rhs,
            SCIPconsIsInitial(cons), SCIPconsIsSeparated(cons), SCIPconsIsEnforced(cons), 
            SCIPconsIsChecked(cons), SCIPconsIsPropagated(cons), SCIPconsIsLocal(cons),
            SCIPconsIsLocal(cons), SCIPconsIsModifiable(cons), 
            SCIPconsIsDynamic(cons), SCIPconsIsRemovable(cons), SCIPconsIsStickingAtNode(cons)) );
   }

   return SCIP_OKAY;
}
#endif




/*
 * constraint specific interface methods
 */

/** creates the handler for exactlp constraints and includes it in SCIP */
SCIP_RETCODE SCIPincludeConshdlrExactlp(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata;

   /* create exactlp constraint handler data */
   SCIP_CALL( conshdlrdataCreate(scip, &conshdlrdata) );

   /* include constraint handler */
   SCIP_CALL( SCIPincludeConshdlr(scip, CONSHDLR_NAME, CONSHDLR_DESC,
         CONSHDLR_SEPAPRIORITY, CONSHDLR_ENFOPRIORITY, CONSHDLR_CHECKPRIORITY,
         CONSHDLR_SEPAFREQ, CONSHDLR_PROPFREQ, CONSHDLR_EAGERFREQ, CONSHDLR_MAXPREROUNDS, 
         CONSHDLR_DELAYSEPA, CONSHDLR_DELAYPROP, CONSHDLR_DELAYPRESOL, CONSHDLR_NEEDSCONS,
         consFreeExactlp, consInitExactlp, consExitExactlp, 
         consInitpreExactlp, consExitpreExactlp, consInitsolExactlp, consExitsolExactlp,
         consDeleteExactlp, consTransExactlp, consInitlpExactlp,
         consSepalpExactlp, consSepasolExactlp, consEnfolpExactlp, consEnfopsExactlp, consCheckExactlp, 
         consPropExactlp, consPresolExactlp, consRespropExactlp, consLockExactlp,
         consActiveExactlp, consDeactiveExactlp, 
         consEnableExactlp, consDisableExactlp,
         consPrintExactlp,
         conshdlrdata) );

#ifdef LINCONSUPGD_PRIORITY
   /* include the linear constraint upgrade in the linear constraint handler */
   SCIP_CALL( SCIPincludeLinconsUpgrade(scip, linconsUpgdExactlp, LINCONSUPGD_PRIORITY) );
#endif

   return SCIP_OKAY;
}

/** creates and captures a exactlp constraint */
SCIP_RETCODE SCIPcreateConsExactlp(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS**           cons,               /**< pointer to hold the created constraint */
   const char*           name,               /**< name of constraint */
   SCIP_OBJSEN           objsense,           /**< objective sense */
   int                   nvars,              /**< number of variables */
   mpq_t*                obj,                /**< objective function values of variables */
   mpq_t*                lb,                 /**< lower bounds of variables */
   mpq_t*                ub,                 /**< upper bounds of variables */
   int                   nconss,              /**< number of constraints */
   mpq_t*                lhs,                /**< left hand sides of constraints */
   mpq_t*                rhs,                /**< right hand sides of constraints */
   int                   nnonz,              /**< number of nonzero elements in the constraint matrix */
   int*                  beg,                /**< start index of each variable in ind- and val-array */
   int*                  ind,                /**< variable indices (var->probindex) of constraint matrix entries */
   mpq_t*                val,                /**< values of nonzero constraint matrix entries */
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
                                              *   are seperated as constraints. */
   SCIP_Bool             removable,          /**< should the relaxation be removed from the LP due to aging or cleanup?
                                              *   Usually set to FALSE. Set to TRUE for 'lazy constraints' and 'user cuts'. */
   SCIP_Bool             stickingatnode      /**< should the constraint always be kept at the node where it was added, even
                                              *   if it may be moved to a more global node?
                                              *   Usually set to FALSE. Set to TRUE to for constraints that represent node data. */
   )
{
   SCIP_CONSHDLR* conshdlr;
   SCIP_CONSDATA* consdata;

   /* find the exactlp constraint handler */
   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);
   if( conshdlr == NULL )
   {
      SCIPerrorMessage("exactlp constraint handler not found\n");
      return SCIP_PLUGINNOTFOUND;
   }

   /* create constraint data */
   SCIP_CALL( consdataCreate(scip, &consdata, objsense, nvars, obj, lb, ub, nconss, lhs, rhs, nnonz, beg, ind, val) );

   /* create constraint */
   SCIP_CALL( SCIPcreateCons(scip, cons, name, conshdlr, consdata, initial, separate, enforce, check, propagate,
         local, modifiable, dynamic, removable, stickingatnode) );

   return SCIP_OKAY;
}
