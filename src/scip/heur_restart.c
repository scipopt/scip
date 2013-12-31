#define SCIP_STATISTIC
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2013 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   heur_restart.c
 * @ingroup PRIMALHEURISTICS
 * @brief  heuristic triggering a restart after addition of cutting planes or variable fixings
 * @author Benjamin Müller
 * @author Ambros Gleixner
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>



#include "scip/heur_restart.h"
#define SCIP_STATISTIC

#define HEUR_NAME             "restart"
#define HEUR_DESC             "heuristic triggering a restart after addition of cutting planes"
#define HEUR_DISPCHAR         '?' /* this character will never be displayed since this is not a primal heuristic */
#define HEUR_PRIORITY         -10000000
#define HEUR_FREQ             -1 /* @TODO set this to -1 */
#define HEUR_FREQOFS          0
#define HEUR_MAXDEPTH         -1
#define HEUR_TIMING           SCIP_HEURTIMING_AFTERNODE
#define HEUR_USESSUBSCIP      FALSE  /**< does the heuristic use a secondary SCIP instance? */

/* default values for user parameters, grouped by parameter type */
#define DEFAULT_MAXRESTARTS   10 /**< maximum number of restarts triggered by this heuristic */

#define DEFAULT_EASYRESTARTS  FALSE /**< use the easy heuristic (as long nruns < nmaxruns => restart!) */

#define DEFAULT_DOMAINRESTARTS  FALSE /**< restart trigger only depends on domain sizes (infinite and finite) */

#define DEFAULT_DOMAINFAC_INFDOMAINS 0.95 /**< trigger a restart if \# inf domains / \# inf domains of the last run < factor */

#define DEFAULT_DOMAINFAC_TOTALDOMAINS 0.95 /**< trigger a restart if total domain / total domain of the last run < factor */

#define DEFAULT_FACTORFIXEDVARS 0.05 /**< if the factor times the number of all variables is less or eq. the number of fixed variables, then apply a restart */

#define DEFAULT_FACTORADDEDCONS 0.10 /**< if the factor times the number of all cons. is less or eq. the number of added cons, then apply a restart */

#define DEFAULT_FACTORSURVCUTS 0.05 /**< if the factor times the number of all cuts is less or eq. the number of survived cuts, then apply a restart */


#define EVENTHDLR_NAME                       "restart_eventhandler"
#define EVENTHDLR_DESC                       "event handler for heur_restart statistics"


#define abs(x) (x > 0) ? x : -x

/*
 * Data structures
 */




/** primal heuristic data */
struct SCIP_HeurData
{
   SCIP_EVENTHDLR*       eventhdlr; /* eventhandler for the statistics */


   /* variables for restart criterions */
   int nmaxruns;                    /* maximum number of restarts */
   SCIP_Bool easyrestart;           /* yes: just restart until we reached nmaxruns   no: more intelligent restart strategy */
   SCIP_Real factorfixedvars;
   SCIP_Real factoraddedcons;
   SCIP_Real factorsurvcuts;

   SCIP_Bool domainrestart;         /* restart trigger only depends on domain sizes (infinite and finite) */
   SCIP_Real factorinfdomains;
   SCIP_Real factortotaldomains;



   SCIP_HASHMAP* addedvarshash;  /* hashmap to identify which variables has been added since the last run */
   SCIP_Real lasttotaldomain;       /* sum of the finite domain sizes */
   int lastninfdomains;             /* total number of variables with infinite domain */


   /* statistic variables */
   int nruns;
   int nfixings;
   int naddedvars;
   int nlastroundfix;
   int ncutstocons;
   int nsurvcuts;

   SCIP_Real avfacinfdomain;
   SCIP_Real avfactotaldomain;
   int nfacinfdomainred;
   int nfactotaldomainred;

   /* statistics for 1. presolve */
   int firstnfixings;
   int firstnvaradds;
   int firstncutstocons;

   SCIP_Real firstfacinfdomain;
   SCIP_Real firstfactotaldomain;


   /* variables to store informations of the last round*/
   int lastnvars;
   int lastncons;
   int lastnaddedcons;

   /* temporary counting variables */
   int tmpnaddedvars;

};

/** event handler data */
struct SCIP_EventhdlrData
{
   SCIP_HEUR* heur;
};



/*
 * Local methods
 */

/**
 * computes the number of pool cuts which will be added to the LP
 */
static
int computeCutsToCons(SCIP* scip)
{
   SCIP_CUT** cuts;
   SCIP_CUT** delayedcuts;
   int ncuts;
   int ndelayedcuts;
   int c;
   int cutstocons;

   /* get all pool cuts */
   ncuts = SCIPgetNPoolCuts(scip);
   cuts = SCIPgetPoolCuts(scip);
   ndelayedcuts = SCIPgetNDelayedPoolCuts(scip);
   delayedcuts = SCIPgetDelayedPoolCuts(scip);


   cutstocons = 0;

   /* check how many pool cuts will be made to linear constraints */
   /*@todo This number is greater than the correct number of pool cuts which will be added! */
   cutstocons = 0;
   for( c = 0; c < ncuts ; ++c )
   {
      SCIP_ROW* row;

      row = SCIPcutGetRow(cuts[c]);
      if( SCIPcutGetAge(cuts[c]) == 0 && SCIProwIsInLP(row) )
         ++cutstocons;
   }

   for( c = 0; c < ndelayedcuts ; ++c )
   {
      SCIP_ROW* row;

      row = SCIPcutGetRow(delayedcuts[c]);
      if( SCIPcutGetAge(delayedcuts[c]) == 0 && SCIProwIsInLP(row) )
         ++cutstocons;
   }

   return cutstocons;

}


/*
 * computes the current sum of finite domain sizes and the number of infinite domains
 * NOTE: the values will be stored in lasttotaldomain and lastninfdomains
 */
static
SCIP_RETCODE computeVariableDomains(SCIP* scip, int* infdomains, SCIP_Real* totaldomain)
{
   SCIP_VAR** vars;
   int nvars;
   int i;

   vars = SCIPgetVars(scip);
   nvars = SCIPgetNVars(scip);

   *infdomains = 0;
   *totaldomain = 0.0;



   for( i = 0; i < nvars; ++i)
   {

      SCIP_Real lb;
      SCIP_Real ub;

      lb = SCIPvarGetLbLocal(vars[i]);
      ub = SCIPvarGetUbLocal(vars[i]);

      if( lb == (-1)*SCIPinfinity(scip)  || ub == SCIPinfinity(scip) ){
         (*infdomains)++;
      }
      else{
         (*totaldomain) += ub - lb;
      }

   }

   return SCIP_OKAY;
}

/*
 * computes the infinite and total domain size of new variables and add them to the hash map
 * NOTE: the values will be stored in lasttotaldomain and lastninfdomains
 */
static
SCIP_RETCODE computeNewVariableDomains(SCIP* scip, SCIP_HASHMAP* addedvarshash, int* infdomains, SCIP_Real* totaldomain)
{
   SCIP_VAR** vars;
   int nvars;
   int i;

   vars = SCIPgetVars(scip);
   nvars = SCIPgetNVars(scip);

   (*infdomains) = 0;
   (*totaldomain) = 0.0;

   for( i = 0; i < nvars; ++i )
   {
      int id;
      id = SCIPvarGetIndex(vars[i]);

      /* true if the variable has been added to the problem */
      if( SCIPhashmapExists(addedvarshash, (void*)(size_t)id ) == FALSE )
      {
         SCIP_Real lb;
         SCIP_Real ub;

         lb = SCIPvarGetLbLocal(vars[i]);
         ub = SCIPvarGetUbLocal(vars[i]);

         if( lb == (-1)*SCIPinfinity(scip)  || ub == SCIPinfinity(scip) )
            (*infdomains)++;

         else
            (*totaldomain) += ub - lb;


         /* now add the variable to the hashmap */
         SCIPhashmapInsert(addedvarshash, (void*)(size_t)id, (void*)(SCIP_Bool)TRUE);

      }
   }

   return SCIP_OKAY;
}


/*
 * Callback methods of primal heuristic
 */

/** execution method of primal heuristic */
static
SCIP_DECL_HEUREXEC(heurExecRestart)
{  /*lint --e{715}*/

   SCIP_HEURDATA* heurdata;
   SCIP_Bool restart;
   int cutstocons;
   SCIP_Real totaldomainall;
   SCIP_Real totaldomainnew;
   int infdomainsall;
   int infdomainsnew;

   SCIP_Real infdomainfac;
   SCIP_Real totaldomainfac;


   /* get heuristic's data */
   heurdata = SCIPheurGetData(heur);

   assert(heurdata != NULL);

   restart = FALSE;


   /*************************************************
    *        criterion for the easy heuristic
    **************************************************/
   if( heurdata->nruns < heurdata->nmaxruns && heurdata->easyrestart == TRUE ){
      restart = TRUE;
   }
   /*************************************************/



   /* ***********************************************
    *     store/compute the statistic informations
    *************************************************/
   cutstocons = computeCutsToCons(scip);

   /* compute domains of all variables and added variables */
   SCIP_CALL( computeVariableDomains(scip, &infdomainsall, &totaldomainall ) );
   SCIP_CALL( computeNewVariableDomains(scip, heurdata->addedvarshash, &infdomainsnew, &totaldomainnew ) );

   infdomainfac = (infdomainsall - infdomainsnew);
   infdomainfac = (heurdata->lastninfdomains == 0) ? 1 : infdomainfac/heurdata->lastninfdomains ;

   totaldomainfac = (heurdata->lasttotaldomain == 0) ? 1 : (totaldomainall - totaldomainnew)/heurdata->lasttotaldomain;



   /* special case for nruns = 0 */
   if( heurdata->nruns == 0 )
   {
      heurdata->firstnfixings = heurdata->lastnvars - SCIPgetNVars(scip) + heurdata->tmpnaddedvars;
      heurdata->firstnvaradds = heurdata->tmpnaddedvars;
      heurdata->firstncutstocons = cutstocons;

      heurdata->firstfacinfdomain = infdomainfac;
      heurdata->firstfactotaldomain = totaldomainfac;
   }





   printf("last inf domain %d \n" , heurdata->lastninfdomains);
   printf("last total domain %e \n" , heurdata->lasttotaldomain);
   printf("INF DOMAIN: NEW %d  OLD %d\n" , infdomainsnew, infdomainsall);
   printf("TOTAL DOMAIN: NEW %e  OLD %e\n" , totaldomainnew, totaldomainall);

   printf("FACTOR inf domain: %e\n" , infdomainfac );
   printf("FACTOR total domain: %e\n" , totaldomainfac );


   /* store the average factors of the domain reduction */
   if( heurdata->nruns > 0 )
   {
      heurdata->avfacinfdomain += infdomainfac;
      heurdata->avfactotaldomain += totaldomainfac;
      heurdata->nfacinfdomainred++;
      heurdata->nfactotaldomainred++;
   }


   heurdata->nfixings += heurdata->lastnvars - SCIPgetNVars(scip) + heurdata->tmpnaddedvars;
   heurdata->naddedvars += heurdata->tmpnaddedvars;
   heurdata->ncutstocons += cutstocons;


   /* how many cuts survived last round? */
   if( SCIPgetNConss(scip) - heurdata->lastncons > 0 )
      heurdata->nsurvcuts += SCIPgetNConss(scip) - heurdata->lastncons;


   /* was there a fixing in this round? */
   if( heurdata->lastnvars - SCIPgetNVars(scip) + heurdata->tmpnaddedvars > 0 )
      heurdata->nlastroundfix = heurdata->nruns;



   /**************************************************************************
    *  triggers for the restarts
    **************************************************************************/
   if( heurdata->easyrestart == FALSE && heurdata->domainrestart == FALSE)
   {
      if( (heurdata->lastnvars - SCIPgetNVars(scip)) >= SCIPgetNVars(scip) * heurdata->factorfixedvars )
         restart = TRUE;
      else if(heurdata->lastnaddedcons >= SCIPgetNConss(scip) *  heurdata->factoraddedcons )
         restart = TRUE;
      else if(  SCIPgetNConss(scip) - heurdata->lastncons >= SCIPgetNConss(scip) *  heurdata->factorsurvcuts)
         restart = TRUE;
      else if( heurdata->nruns == 0 )
         restart = TRUE;
   }

   if( heurdata->domainrestart == TRUE )
   {
      if( infdomainfac <= heurdata->factorinfdomains )
         restart = TRUE;

      if( totaldomainfac <= heurdata->factortotaldomains )
         restart = TRUE;

      if( heurdata->nruns == 0 )
         restart = TRUE;
   }


   /*****************************************************************************/



   /* update some variables for the next run */
   heurdata->tmpnaddedvars = 0;
   heurdata->lastnvars = SCIPgetNVars(scip);
   heurdata->lastncons = SCIPgetNConss(scip);
   heurdata->lastnaddedcons = cutstocons;

   heurdata->lastninfdomains = infdomainsall;
   heurdata->lasttotaldomain = totaldomainall;

   /* do the restart if necessary */
   if( restart == TRUE )
   {
      SCIP_CALL( SCIPrestartSolve(scip) );
      heurdata->nruns++;
   }

   return SCIP_OKAY;
}


/** deinitialization method of primal heuristic (called before transformed problem is freed) */
static
SCIP_DECL_HEUREXIT(heurExitRestart)
{  /*lint --e{715}*/

   SCIP_HEURDATA* heurdata;
   heurdata = SCIPheurGetData(heur);

   assert(heurdata != NULL);



   /* print the statistics */
   SCIPstatisticMessage("##NRuns: %d | NAddedCuts: %d | NSurvCuts: %d | NFixVars: %d | LastFixRound: %d | NAddedVars: %d | FirstFixings: %d | FirstVarAdds: %d | FirstAddedCuts: %d | FirstFacInfDom: %e | FirstFacTotalDom: %e | AvFacInfDom: %e | num :%d | AvFacTotalDom: %e | num: %d |\n", \
      heurdata->nruns, heurdata->ncutstocons, heurdata->nsurvcuts, heurdata->nfixings, heurdata->nlastroundfix, heurdata->naddedvars, heurdata->firstnfixings, heurdata->firstnvaradds, heurdata->firstncutstocons, \
      heurdata->firstfacinfdomain, heurdata->firstfactotaldomain, heurdata->avfacinfdomain/heurdata->nfacinfdomainred , heurdata->nfacinfdomainred, heurdata->avfactotaldomain/heurdata->nfactotaldomainred, heurdata->nfactotaldomainred );


   /* free the hash map */
   SCIPhashmapFree(&heurdata->addedvarshash);

   /* now free heurdata */
   SCIPfreeMemory(scip, &heurdata);
   SCIPheurSetData(heur, NULL);

   return SCIP_OKAY;
}


/** initialization method of primal heuristic (called after problem was transformed) */
static
SCIP_DECL_HEURINIT(heurInitRestart)
{  /*lint --e{715}*/
   SCIP_HEURDATA* heurdata;

   /* get heuristic's data */
   heurdata = SCIPheurGetData(heur);

   assert(heurdata != NULL);

   /* alloc memory for the hashtable */
   SCIP_CALL( SCIPhashmapCreate(&heurdata->addedvarshash, SCIPblkmem(scip), 10*SCIPgetNVars(scip) ) );

   /* initialize heurdata variables */
   heurdata->lastnvars = SCIPgetNVars(scip);
   heurdata->lastncons = SCIPgetNConss(scip);
   heurdata->nlastroundfix = -1;
   heurdata->avfacinfdomain = 0;
   heurdata->avfactotaldomain = 0;
   heurdata->nfacinfdomainred = 0;
   heurdata->nfactotaldomainred = 0;

   /* when a variable has been added */
   SCIP_CALL( SCIPcatchEvent( scip, SCIP_EVENTTYPE_VARADDED, heurdata->eventhdlr, NULL, NULL) );

   /* initialize the domain variables for the statistics */
   heurdata->lastninfdomains = 0;
   heurdata->lasttotaldomain = 0.0;
   SCIP_CALL( computeNewVariableDomains(scip, heurdata->addedvarshash, &heurdata->lastninfdomains, &heurdata->lasttotaldomain) );

   printf("inf domains: %d total domain size: %e\n", heurdata->lastninfdomains, heurdata->lasttotaldomain );

   return SCIP_OKAY;
}



/*
 * Callback methods of event handler
 */

static
SCIP_DECL_EVENTEXEC(eventExecRestartStatistics)
{  /*lint --e{715}*/

   SCIP_HEURDATA* heurdata;
   SCIP_EVENTHDLRDATA* eventhdlrdata;

   assert(scip != NULL);

   eventhdlrdata = SCIPeventhdlrGetData(eventhdlr);
   heurdata = SCIPheurGetData(eventhdlrdata->heur);

   assert(eventhdlrdata != NULL);
   assert(heurdata != NULL);

   /* count the number of added variables */
   if( SCIPeventGetType(event) ==  SCIP_EVENTTYPE_VARADDED )
   {
      heurdata->tmpnaddedvars++;
   }
   return SCIP_OKAY;
}

static
SCIP_DECL_EVENTFREE(eventFree)
{
   SCIP_EVENTHDLRDATA* eventhdlrdata;
   eventhdlrdata = SCIPeventhdlrGetData(eventhdlr);
   assert(eventhdlrdata != NULL);
   SCIPfreeMemory(scip, &eventhdlrdata);
   SCIPeventhdlrSetData(eventhdlr, NULL);
   return SCIP_OKAY;
}




/*
 * primal heuristic specific interface methods
 */

/** creates the restart primal heuristic and includes it in SCIP */
SCIP_RETCODE SCIPincludeHeurRestart(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_HEURDATA* heurdata;
   SCIP_HEUR* heur;
   SCIP_EVENTHDLRDATA* eventhdlrdata;

   /* create restart primal heuristic data */
   SCIP_CALL( SCIPallocMemory(scip, &heurdata) );

   /* create eventhandler data */
   SCIP_CALL( SCIPallocMemory(scip, &eventhdlrdata) );


   /* include primal heuristic */
   SCIP_CALL( SCIPincludeHeurBasic(scip, &heur,
         HEUR_NAME, HEUR_DESC, HEUR_DISPCHAR, HEUR_PRIORITY, HEUR_FREQ, HEUR_FREQOFS,
         HEUR_MAXDEPTH, HEUR_TIMING, HEUR_USESSUBSCIP, heurExecRestart, heurdata) );

   /* set non fundamental callbacks via setter functions */
   SCIP_CALL( SCIPsetHeurExit(scip, heur, heurExitRestart) );
   SCIP_CALL( SCIPsetHeurInit(scip, heur, heurInitRestart) );

   assert(heur != NULL);

   /* add restart primal heuristic parameters */
   SCIP_CALL( SCIPaddIntParam(scip, "heuristics/"HEUR_NAME"/maxrestarts",
         "maximum number of restarts triggered by this heuristic",
         &heurdata->nmaxruns, FALSE, DEFAULT_MAXRESTARTS, 0, INT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "heuristics/"HEUR_NAME"/easyrestarts",
         "yes: just restart until we reached nmaxruns   no: more intelligent restart strategy",
         &heurdata->easyrestart, FALSE, DEFAULT_EASYRESTARTS, 0, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "heuristics/"HEUR_NAME"/domainrestarts",
         "restart trigger only depends on domain sizes (infinite and finite)",
         &heurdata->domainrestart, FALSE, DEFAULT_DOMAINRESTARTS, 0, NULL) );


   SCIP_CALL( SCIPaddRealParam(scip, "heuristics/"HEUR_NAME"/factorinfdomains",
         "trigger a restart if #inf domains / # inf domains of the last run < factor",
         &heurdata->factorinfdomains, FALSE, DEFAULT_DOMAINFAC_INFDOMAINS, 0, 1, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip, "heuristics/"HEUR_NAME"/factortotaldomains",
         "trigger a restart if total domain / total domain of the last run < factor",
         &heurdata->factortotaldomains, FALSE, DEFAULT_DOMAINFAC_TOTALDOMAINS, 0, 1, NULL, NULL) );


   SCIP_CALL( SCIPaddRealParam(scip, "heuristics/"HEUR_NAME"/factorfixedvars",
         "if the factor times the number of all variables is less or eq. the number of fixed variables, then apply a restart",
         &heurdata->factorfixedvars, FALSE, DEFAULT_FACTORFIXEDVARS, 0, 1, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip, "heuristics/"HEUR_NAME"/factoraddedcons",
         "if the factor times the number of all cons. is less or eq. the number of added cons, then apply a restart",
         &heurdata->factoraddedcons, FALSE, DEFAULT_FACTORADDEDCONS, 0, 1, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip, "heuristics/"HEUR_NAME"/factorsurvcuts",
         " if the factor times the number of all cuts is less or eq. the number of survived cuts, then apply a restart",
         &heurdata->factorsurvcuts, FALSE, DEFAULT_FACTORSURVCUTS, 0, 1, NULL, NULL) );




   /* initialize eventhdlr data */
   eventhdlrdata->heur = heur;

   /* include event handler */
   SCIP_CALL( SCIPincludeEventhdlrBasic(scip, &heurdata->eventhdlr, EVENTHDLR_NAME, EVENTHDLR_DESC, eventExecRestartStatistics, eventhdlrdata) );
   SCIP_CALL( SCIPsetEventhdlrFree(scip, heurdata->eventhdlr, eventFree) );

   return SCIP_OKAY;
}
