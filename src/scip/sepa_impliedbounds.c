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
#pragma ident "@(#) $Id: sepa_impliedbounds.c,v 1.1 2005/04/11 10:56:15 bzfpfend Exp $"

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
   Real*            varsolvals,         /**< LP solution of active problem variables */
   Real             val1,               /**< given coefficient of first variable */
   VAR*             var1,               /**< given first variable */
   Real             val2,               /**< given coefficient of second variable */
   VAR*             var2,               /**< given second variable */
   Real             rhs,                /**< given right hand side of the cut to add */
   int*             ncuts               /**< pointer to update number of cuts added */
   )
{
   Real activity;

   assert(scip != NULL);

   activity = val1 * varsolvals[SCIPvarGetProbindex(var1)] + val2 * varsolvals[SCIPvarGetProbindex(var2)]; 
   debugMessage(" -> %g<%s> + %g<%s> <= %g (act: %g)\n", val1, SCIPvarGetName(var1), val2, SCIPvarGetName(var2),
      rhs, activity);

   if( SCIPisEfficacious(scip, activity - rhs) )
   {
      ROW* cut;
      char cutname[MAXSTRLEN];
      
      /* creates the cut */
      sprintf(cutname, "implbd%d_%d", SCIPgetNLPs(scip), *ncuts);
      CHECK_OKAY( SCIPcreateEmptyRow(scip, &cut, cutname, -SCIPinfinity(scip), rhs, FALSE, FALSE, TRUE) );
      
      CHECK_OKAY( SCIPcacheRowExtensions(scip, cut) );
 
      CHECK_OKAY( SCIPaddVarToRow(scip, cut, var1, val1) );
      CHECK_OKAY( SCIPaddVarToRow(scip, cut, var2, val2) );
      
      CHECK_OKAY( SCIPflushRowExtensions(scip, cut) );
      
#ifdef DEBUG
      printf(" -> found cut (activity = %g): ", activity);
      SCIPprintRow(scip, cut, NULL);
#endif
      /* adds the cut */
      CHECK_OKAY( SCIPaddCut(scip, cut, FALSE) );
      (*ncuts)++;
      
      /* release the row */
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
   Real* varsolvals;
   int nbinvars;
   int nvars;
   int ncuts;
   int i;
   int j;
   BOUNDTYPE* impltypes; 
   Real* implbounds; 
   VAR** implvars;
   int nimpl;

   assert(sepa != NULL);
   assert(scip != NULL);
 
   *result = SCIP_DIDNOTRUN;

   /* gets active problem variables */
   CHECK_OKAY( SCIPgetVarsData(scip, &vars, &nvars, &nbinvars, NULL, NULL, NULL) );
   assert(vars != NULL);

   /* there is nothing to do, if there are no binary variables */
   if( nbinvars == 0 )
      return SCIP_OKAY;

   *result = SCIP_DIDNOTFIND;
   ncuts = 0;

   debugMessage("searching for implied bound cuts\n");

   /* gets data structure */
   CHECK_OKAY( SCIPallocBufferArray(scip, &varsolvals, nvars) );

   /* gets LP solution of all active variables */
   CHECK_OKAY( SCIPgetVarSols(scip, nvars, vars, varsolvals) );

   /* search binary variables for violated implications */
   for( i = 0; i < nbinvars; i++ )
   {
      /* gets all implications x >= 1 ==> y <= p or y >= p */
      nimpl = SCIPvarGetNBinImpls(vars[i], TRUE);
      implvars = SCIPvarGetImplVars(vars[i], TRUE);
      impltypes = SCIPvarGetImplTypes(vars[i], TRUE);
      implbounds = SCIPvarGetImplBounds(vars[i], TRUE);
      
      /* tries to add cuts for implications of x >= 1  
       *    y <= ub + x * (p - ub)  <==>  y + (ub - p) * x <=  ub          if x >= 1 ==> y <= p,   
       *    y >= lb + x * (p - lb)  <==> -y + (p - lb) * x <= -lb          if x >= 1 ==> y >= p
       * with lb (ub) global lower (upper) bound of y
       */
      for( j = 0; j < nimpl; j++ )
      {
         assert(SCIPvarGetType(vars[i]) == SCIP_VARTYPE_BINARY);

         if( impltypes[j] == SCIP_BOUNDTYPE_UPPER )
         {
            Real ub;
         
            /* implication x >= 1 ==> y <= p */
            ub = SCIPvarGetUbGlobal(implvars[j]);
            assert(SCIPisLE(scip, implbounds[j], ub));
            
            /* adds cut if violated */
            CHECK_OKAY( addCut(scip, varsolvals, 1.0, implvars[j], (ub - implbounds[j]), vars[i], ub, &ncuts) );
         }
         else
         {
            Real lb;

            /* implication x >= 1 ==> y >= p */
            lb = SCIPvarGetLbGlobal(implvars[j]);
            assert(impltypes[j] == SCIP_BOUNDTYPE_LOWER && SCIPisGE(scip, implbounds[j], lb));
 
            /* adds cut if violated */
            CHECK_OKAY( addCut(scip, varsolvals, -1.0, implvars[j], (implbounds[j] - lb), vars[i], -lb, &ncuts) );
         }
      } 

      /* gets all implications x <= 0 ==> y <= p or y >= p */
      nimpl = SCIPvarGetNBinImpls(vars[i], FALSE);
      implvars = SCIPvarGetImplVars(vars[i], FALSE);
      impltypes = SCIPvarGetImplTypes(vars[i], FALSE);
      implbounds = SCIPvarGetImplBounds(vars[i], FALSE);
      
      /* tries to add cuts for implications of x >= 1  
       *    y <= p + x * (ub - p)  <==>  y + (p - ub) * x <=  p            if x <= 0 ==> y <= p,   
       *    y >= p + x * (lb - p)  <==> -y + (lb - p) * x <= -p            if x <= 0 ==> y >= p
       * with lb (ub) global lower (upper) bound of y
       */
      for( j = 0; j < nimpl; j++ )
      {
         assert(SCIPvarGetType(vars[i]) == SCIP_VARTYPE_BINARY);

         if( impltypes[j] == SCIP_BOUNDTYPE_UPPER )
         {
            Real ub;
         
            /* implication x <= 0 ==> y <= p */
            ub = SCIPvarGetUbGlobal(implvars[j]);
            assert(SCIPisLE(scip, implbounds[j], ub));
 
            /* adds cut if violated */
            CHECK_OKAY( addCut(scip, varsolvals, 1.0, implvars[j], (implbounds[j] - ub), vars[i], implbounds[j], &ncuts) );
         }
         else
         {
            Real lb;

            /* implication x <= 0 ==> y >= p */
            lb = SCIPvarGetLbGlobal(implvars[j]);
            assert(impltypes[j] == SCIP_BOUNDTYPE_LOWER && SCIPisGE(scip, implbounds[j], lb));
 
            /* adds cut if violated */
            CHECK_OKAY( addCut(scip, varsolvals, -1.0, implvars[j], (lb - implbounds[j]), vars[i], -implbounds[j],
                  &ncuts) );
         }
      } 
   }

   if( ncuts > 0 )
      *result = SCIP_SEPARATED;

   /* frees data structures */
   SCIPfreeBufferArray(scip, &varsolvals);

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
