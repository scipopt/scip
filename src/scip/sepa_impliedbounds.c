/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2005 Tobias Achterberg                              */
/*                                                                           */
/*                  2002-2005 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the SCIP Academic License.        */
/*                                                                           */
/*  You should have received a copy of the SCIP Academic License             */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#pragma ident "@(#) $Id: sepa_impliedbounds.c,v 1.2 2005/04/12 16:56:17 bzfpfend Exp $"

/**@file   sepa_impliedbounds.c
 * @brief  implied bounds separator
 * @author Kati Wolter
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>

#include "scip/sepa_impliedbounds.h"


#define SEPA_NAME              "impliedbounds"
#define SEPA_DESC              "implied bounds separator"
#define SEPA_PRIORITY               -50
#define SEPA_FREQ                    10
#define SEPA_DELAY                FALSE /**< should separation method be delayed, if other separators found cuts? */




/*
 * Data structures
 */

/** separator data */
struct SepaData
{
};




/*
 * Local methods
 */

/* adds given cut with two variables, if it is violated */
static 
RETCODE addCut(
   SCIP*            scip,               /**< SCIP data structure */
   Real             val1,               /**< given coefficient of first variable */
   VAR*             var1,               /**< given first variable */
   Real             solval1,            /**< current LP solution value of first variable */
   Real             val2,               /**< given coefficient of second variable */
   VAR*             var2,               /**< given second variable */
   Real             solval2,            /**< current LP solution value of second variable */
   Real             rhs,                /**< given right hand side of the cut to add */
   int*             ncuts               /**< pointer to update number of cuts added */
   )
{
   Real activity;

   assert(ncuts != NULL);

   /* calculate activity of cut */
   activity = val1 * solval1 + val2 * solval2;
   /*debugMessage(" -> %g<%s>[%g] + %g<%s>[%g] <= %g (act: %g)\n", 
     val1, SCIPvarGetName(var1), solval1, val2, SCIPvarGetName(var2), solval2, rhs, activity);*/

   /* check, if cut is violated */
   if( SCIPisEfficacious(scip, activity - rhs) )
   {
      ROW* cut;
      char cutname[MAXSTRLEN];

      /* create cut */
      sprintf(cutname, "implbd%d_%d", SCIPgetNLPs(scip), *ncuts);
      CHECK_OKAY( SCIPcreateEmptyRow(scip, &cut, cutname, -SCIPinfinity(scip), rhs, FALSE, FALSE, TRUE) );
      CHECK_OKAY( SCIPcacheRowExtensions(scip, cut) );
      CHECK_OKAY( SCIPaddVarToRow(scip, cut, var1, val1) );
      CHECK_OKAY( SCIPaddVarToRow(scip, cut, var2, val2) );
      CHECK_OKAY( SCIPflushRowExtensions(scip, cut) );
      
#ifdef DEBUG
      debugMessage(" -> found cut (activity = %g): ", activity);
      SCIPprintRow(scip, cut, NULL);
#endif

      /* add cut */
      CHECK_OKAY( SCIPaddCut(scip, cut, FALSE) );
      (*ncuts)++;
      
      /* release cut */
      CHECK_OKAY( SCIPreleaseRow(scip, &cut) );
   }
   
   return SCIP_OKAY;
}


/*
 * Callback methods of separator
 */


/** destructor of separator to free user data (called when SCIP is exiting) */
#define sepaFreeImpliedbounds NULL



/** initialization method of separator (called after problem was transformed) */
#define sepaInitImpliedbounds NULL


/** deinitialization method of separator (called before transformed problem is freed) */
#define sepaExitImpliedbounds NULL


/** solving process initialization method of separator (called when branch and bound process is about to begin) */
#define sepaInitsolImpliedbounds NULL


/** solving process deinitialization method of separator (called before branch and bound process data is freed) */
#define sepaExitsolImpliedbounds NULL


/** execution method of separator */
static
DECL_SEPAEXEC(sepaExecImpliedbounds)
{  /*lint --e{715}*/
   VAR** vars;
   VAR** fracvars;
   Real* solvals;
   Real* fracvals;
   int nvars;
   int nbinvars;
   int nfracs;
   int ncuts;
   int i;

   assert(sepa != NULL);
   assert(scip != NULL);
 
   *result = SCIP_DIDNOTRUN;

   /* gets active problem variables */
   CHECK_OKAY( SCIPgetVarsData(scip, &vars, &nvars, &nbinvars, NULL, NULL, NULL) );
   if( nbinvars == 0 )
      return SCIP_OKAY;

   /* get fractional problem variables */
   CHECK_OKAY( SCIPgetLPBranchCands(scip, &fracvars, &fracvals, NULL, &nfracs, NULL) );
   if( nfracs == 0 )
      return SCIP_OKAY;

   /* get solution values for all variables */
   CHECK_OKAY( SCIPallocBufferArray(scip, &solvals, nvars) );
   CHECK_OKAY( SCIPgetVarSols(scip, nvars, vars, solvals) );

   *result = SCIP_DIDNOTFIND;
   ncuts = 0;

   debugMessage("searching for implied bound cuts\n");

   /* search binary variables for violated implications */
   for( i = 0; i < nfracs; i++ )
   {
      BOUNDTYPE* impltypes; 
      Real* implbounds; 
      VAR** implvars;
      int nimpl;
      int j;

      /* only process binary variables */
      if( SCIPvarGetType(fracvars[i]) != SCIP_VARTYPE_BINARY )
         continue;

      /* get implications of x == 1 */
      nimpl = SCIPvarGetNImpls(fracvars[i], TRUE);
      implvars = SCIPvarGetImplVars(fracvars[i], TRUE);
      impltypes = SCIPvarGetImplTypes(fracvars[i], TRUE);
      implbounds = SCIPvarGetImplBounds(fracvars[i], TRUE);
      
      /*debugMessage("%d implications for <%s>[%g] == 1\n", nimpl, SCIPvarGetName(fracvars[i]), fracvals[i]);*/

      /* try to add cuts for implications of x == 1
       *    x == 1 -> y <= p:  y <= ub + x * (p - ub)  <==>  y + (ub - p) * x <=  ub
       *    x == 1 -> y >= p:  y >= lb + x * (p - lb)  <==> -y + (p - lb) * x <= -lb
       * with lb (ub) global lower (upper) bound of y
       */
      for( j = 0; j < nimpl; j++ )
      {
         Real solval;

         solval = solvals[SCIPvarGetProbindex(implvars[j])];
         if( impltypes[j] == SCIP_BOUNDTYPE_UPPER )
         {
            Real ub;
         
            /* implication x == 1 -> y <= p */
            ub = SCIPvarGetUbGlobal(implvars[j]);
            assert(SCIPisLE(scip, implbounds[j], ub));
            
            /* add cut if violated */
            CHECK_OKAY( addCut(scip, 1.0, implvars[j], solval, (ub - implbounds[j]), fracvars[i], fracvals[i],
                  ub, &ncuts) );
         }
         else
         {
            Real lb;

            /* implication x == 1 -> y >= p */
            lb = SCIPvarGetLbGlobal(implvars[j]);
            assert(impltypes[j] == SCIP_BOUNDTYPE_LOWER && SCIPisGE(scip, implbounds[j], lb));
 
            /* add cut if violated */
            CHECK_OKAY( addCut(scip, -1.0, implvars[j], solval, (implbounds[j] - lb), fracvars[i], fracvals[i],
                  -lb, &ncuts) );
         }
      } 

      /* get implications of x == 0 */
      nimpl = SCIPvarGetNImpls(fracvars[i], FALSE);
      implvars = SCIPvarGetImplVars(fracvars[i], FALSE);
      impltypes = SCIPvarGetImplTypes(fracvars[i], FALSE);
      implbounds = SCIPvarGetImplBounds(fracvars[i], FALSE);
      
      /*debugMessage("%d implications for <%s>[%g] == 0\n", nimpl, SCIPvarGetName(fracvars[i]), fracvals[i]);*/

      /* try to add cuts for implications of x == 0
       *    x == 0 -> y <= p:  y <= p + x * (ub - p)  <==>  y + (p - ub) * x <=  p
       *    x == 0 -> y >= p:  y >= p + x * (lb - p)  <==> -y + (lb - p) * x <= -p
       * with lb (ub) global lower (upper) bound of y
       */
      for( j = 0; j < nimpl; j++ )
      {
         Real solval;

         solval = solvals[SCIPvarGetProbindex(implvars[j])];
         if( impltypes[j] == SCIP_BOUNDTYPE_UPPER )
         {
            Real ub;
         
            /* implication x == 0 -> y <= p */
            ub = SCIPvarGetUbGlobal(implvars[j]);
            assert(SCIPisLE(scip, implbounds[j], ub));
 
            /* add cut if violated */
            CHECK_OKAY( addCut(scip, 1.0, implvars[j], solval, (implbounds[j] - ub), fracvars[i], fracvals[i],
                  implbounds[j], &ncuts) );
         }
         else
         {
            Real lb;

            /* implication x == 0 -> y >= p */
            lb = SCIPvarGetLbGlobal(implvars[j]);
            assert(impltypes[j] == SCIP_BOUNDTYPE_LOWER && SCIPisGE(scip, implbounds[j], lb));
 
            /* add cut if violated */
            CHECK_OKAY( addCut(scip, -1.0, implvars[j], solval, (lb - implbounds[j]), fracvars[i], fracvals[i],
                  -implbounds[j], &ncuts) );
         }
      }
   }

   /* adjust result code */
   if( ncuts > 0 )
      *result = SCIP_SEPARATED;

   /* free temporary memory */
   SCIPfreeBufferArray(scip, &solvals);

   return SCIP_OKAY;
}





/*
 * separator specific interface methods
 */

/** creates the impliedbounds separator and includes it in SCIP */
RETCODE SCIPincludeSepaImpliedbounds(
   SCIP*            scip                /**< SCIP data structure */
   )
{
   SEPADATA* sepadata;

   /* create impliedbounds separator data */
   sepadata = NULL;

   /* include separator */
   CHECK_OKAY( SCIPincludeSepa(scip, SEPA_NAME, SEPA_DESC, SEPA_PRIORITY, SEPA_FREQ, SEPA_DELAY,
         sepaFreeImpliedbounds, sepaInitImpliedbounds, sepaExitImpliedbounds, 
         sepaInitsolImpliedbounds, sepaExitsolImpliedbounds, sepaExecImpliedbounds,
         sepadata) );

   return SCIP_OKAY;
}
