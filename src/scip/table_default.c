/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2017 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   table_default.c
 * @brief  default statistics tables
 * @author Tristan Gally
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>

#include "scip/table_default.h"

#define TABLE_NAME_STATUS                "status"
#define TABLE_DESC_STATUS                "status statistics table"
#define TABLE_POSITION_STATUS            0                      /**< the position of the statistics table */
#define TABLE_EARLIEST_STAGE_STATUS      SCIP_STAGE_INIT        /**< output of the statistics table is only printed from this stage onwards */

#define TABLE_NAME_TIMING                "timing"
#define TABLE_DESC_TIMING                "timing statistics table"
#define TABLE_POSITION_TIMING            1000                   /**< the position of the statistics table */
#define TABLE_EARLIEST_STAGE_TIMING      SCIP_STAGE_PROBLEM     /**< output of the statistics table is only printed from this stage onwards */

#define TABLE_NAME_ORIGPROB              "original problem"
#define TABLE_DESC_ORIGPROB              "original problem statistics table"
#define TABLE_POSITION_ORIGPROB          2000                   /**< the position of the statistics table */
#define TABLE_EARLIEST_STAGE_ORIGPROB    SCIP_STAGE_PROBLEM     /**< output of the statistics table is only printed from this stage onwards */

#define TABLE_NAME_TRANSPROB             "presolved problem"
#define TABLE_DESC_TRANSPROB             "presolved problem statistics table"
#define TABLE_POSITION_TRANSPROB         3000                   /**< the position of the statistics table */
#define TABLE_EARLIEST_STAGE_TRANSPROB   SCIP_STAGE_TRANSFORMED /**< output of the statistics table is only printed from this stage onwards */

#define TABLE_NAME_PRESOL                "presolver"
#define TABLE_DESC_PRESOL                "presolver statistics table"
#define TABLE_POSITION_PRESOL            4000                   /**< the position of the statistics table */
#define TABLE_EARLIEST_STAGE_PRESOL      SCIP_STAGE_TRANSFORMED /**< output of the statistics table is only printed from this stage onwards */

#define TABLE_NAME_CONS                  "constraint"
#define TABLE_DESC_CONS                  "constraint statistics table"
#define TABLE_POSITION_CONS              5000                   /**< the position of the statistics table */
#define TABLE_EARLIEST_STAGE_CONS        SCIP_STAGE_TRANSFORMED /**< output of the statistics table is only printed from this stage onwards */

#define TABLE_NAME_CONSTIMING            "constraint timing"
#define TABLE_DESC_CONSTIMING            "constraint timing statistics table"
#define TABLE_POSITION_CONSTIMING        6000                   /**< the position of the statistics table */
#define TABLE_EARLIEST_STAGE_CONSTIMING  SCIP_STAGE_TRANSFORMED /**< output of the statistics table is only printed from this stage onwards */

#define TABLE_NAME_PROP                  "propagator"
#define TABLE_DESC_PROP                  "propagator statistics table"
#define TABLE_POSITION_PROP              7000                   /**< the position of the statistics table */
#define TABLE_EARLIEST_STAGE_PROP        SCIP_STAGE_TRANSFORMED /**< output of the statistics table is only printed from this stage onwards */

#define TABLE_NAME_CONFLICT              "conflict"
#define TABLE_DESC_CONFLICT              "conflict statistics table"
#define TABLE_POSITION_CONFLICT          8000                   /**< the position of the statistics table */
#define TABLE_EARLIEST_STAGE_CONFLICT    SCIP_STAGE_TRANSFORMED /**< output of the statistics table is only printed from this stage onwards */

#define TABLE_NAME_SEPA                  "separator"
#define TABLE_DESC_SEPA                  "separator statistics table"
#define TABLE_POSITION_SEPA              9000                   /**< the position of the statistics table */
#define TABLE_EARLIEST_STAGE_SEPA        SCIP_STAGE_SOLVING     /**< output of the statistics table is only printed from this stage onwards */

#define TABLE_NAME_PRICER                "pricer"
#define TABLE_DESC_PRICER                "pricer statistics table"
#define TABLE_POSITION_PRICER            10000                  /**< the position of the statistics table */
#define TABLE_EARLIEST_STAGE_PRICER      SCIP_STAGE_SOLVING     /**< output of the statistics table is only printed from this stage onwards */

#define TABLE_NAME_BRANCH                "branching rules"
#define TABLE_DESC_BRANCH                "branching rules statistics table"
#define TABLE_POSITION_BRANCH            11000                  /**< the position of the statistics table */
#define TABLE_EARLIEST_STAGE_BRANCH      SCIP_STAGE_SOLVING     /**< output of the statistics table is only printed from this stage onwards */

#define TABLE_NAME_HEUR                  "heuristics"
#define TABLE_DESC_HEUR                  "heuristics statistics table"
#define TABLE_POSITION_HEUR              12000                  /**< the position of the statistics table */
#define TABLE_EARLIEST_STAGE_HEUR        SCIP_STAGE_PRESOLVING  /**< output of the statistics table is only printed from this stage onwards */

#define TABLE_NAME_COMPRESSION           "compression"
#define TABLE_DESC_COMPRESSION           "compression statistics table"
#define TABLE_POSITION_COMPRESSION       13000                  /**< the position of the statistics table */
#define TABLE_EARLIEST_STAGE_COMPRESSION SCIP_STAGE_PRESOLVING  /**< output of the statistics table is only printed from this stage onwards */

#define TABLE_NAME_LP                    "lp"
#define TABLE_DESC_LP                    "lp statistics table"
#define TABLE_POSITION_LP                14000                  /**< the position of the statistics table */
#define TABLE_EARLIEST_STAGE_LP          SCIP_STAGE_SOLVING     /**< output of the statistics table is only printed from this stage onwards */

#define TABLE_NAME_NLP                   "nlp"
#define TABLE_DESC_NLP                   "nlp statistics table"
#define TABLE_POSITION_NLP               15000                  /**< the position of the statistics table */
#define TABLE_EARLIEST_STAGE_NLP         SCIP_STAGE_SOLVING     /**< output of the statistics table is only printed from this stage onwards */

#define TABLE_NAME_RELAX                 "relaxator"
#define TABLE_DESC_RELAX                 "relaxator statistics table"
#define TABLE_POSITION_RELAX             16000                  /**< the position of the statistics table */
#define TABLE_EARLIEST_STAGE_RELAX       SCIP_STAGE_SOLVING     /**< output of the statistics table is only printed from this stage onwards */

#define TABLE_NAME_TREE                  "tree"
#define TABLE_DESC_TREE                  "tree statistics table"
#define TABLE_POSITION_TREE              17000                  /**< the position of the statistics table */
#define TABLE_EARLIEST_STAGE_TREE        SCIP_STAGE_SOLVING     /**< output of the statistics table is only printed from this stage onwards */

#define TABLE_NAME_ROOT                  "root"
#define TABLE_DESC_ROOT                  "root statistics table"
#define TABLE_POSITION_ROOT              18000                  /**< the position of the statistics table */
#define TABLE_EARLIEST_STAGE_ROOT        SCIP_STAGE_SOLVING     /**< output of the statistics table is only printed from this stage onwards */

#define TABLE_NAME_SOL                   "solution"
#define TABLE_DESC_SOL                   "solution statistics table"
#define TABLE_POSITION_SOL               19000                  /**< the position of the statistics table */
#define TABLE_EARLIEST_STAGE_SOL         SCIP_STAGE_PRESOLVING  /**< output of the statistics table is only printed from this stage onwards */

#define TABLE_NAME_CONC                  "concurrent solver"
#define TABLE_DESC_CONC                  "concurrent solver statistics table"
#define TABLE_POSITION_CONC              20000                  /**< the position of the statistics table */
#define TABLE_EARLIEST_STAGE_CONC        SCIP_STAGE_TRANSFORMED /**< output of the statistics table is only printed from this stage onwards */

/*
 * Callback methods of statistics table
 */

/** copy method for statistics table plugins (called when SCIP copies plugins) */
static
SCIP_DECL_TABLECOPY(tableCopyDefault)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(table != NULL);

   /* call inclusion method of statistics tables */
   SCIP_CALL( SCIPincludeTableDefault(scip) );

   return SCIP_OKAY;
}

/** output method of statistics table to output file stream 'file' */
static
SCIP_DECL_TABLEOUTPUT(tableOutputStatus)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(table != NULL);

   SCIPprintStatusStatistics(scip, file);

   return SCIP_OKAY;
}

/** output method of statistics table to output file stream 'file' */
static
SCIP_DECL_TABLEOUTPUT(tableOutputTiming)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(table != NULL);

   SCIPprintTimingStatistics(scip, file);

   return SCIP_OKAY;
}

/** output method of statistics table to output file stream 'file' */
static
SCIP_DECL_TABLEOUTPUT(tableOutputOrigProb)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(table != NULL);

   SCIPprintOrigProblemStatistics(scip, file);

   return SCIP_OKAY;
}

/** output method of statistics table to output file stream 'file' */
static
SCIP_DECL_TABLEOUTPUT(tableOutputTransProb)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(table != NULL);

   SCIPprintTransProblemStatistics(scip, file);

   return SCIP_OKAY;
}

/** output method of statistics table to output file stream 'file' */
static
SCIP_DECL_TABLEOUTPUT(tableOutputPresol)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(table != NULL);

   SCIPprintPresolverStatistics(scip, file);

   return SCIP_OKAY;
}

/** output method of statistics table to output file stream 'file' */
static
SCIP_DECL_TABLEOUTPUT(tableOutputCons)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(table != NULL);

   SCIPprintConstraintStatistics(scip, file);

   return SCIP_OKAY;
}

/** output method of statistics table to output file stream 'file' */
static
SCIP_DECL_TABLEOUTPUT(tableOutputConstiming)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(table != NULL);

   SCIPprintConstraintTimingStatistics(scip, file);

   return SCIP_OKAY;
}

/** output method of statistics table to output file stream 'file' */
static
SCIP_DECL_TABLEOUTPUT(tableOutputProp)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(table != NULL);

   SCIPprintPropagatorStatistics(scip, file);

   return SCIP_OKAY;
}

/** output method of statistics table to output file stream 'file' */
static
SCIP_DECL_TABLEOUTPUT(tableOutputConflict)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(table != NULL);

   SCIPprintConflictStatistics(scip, file);

   return SCIP_OKAY;
}

/** output method of statistics table to output file stream 'file' */
static
SCIP_DECL_TABLEOUTPUT(tableOutputSepa)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(table != NULL);

   SCIPprintSeparatorStatistics(scip, file);

   return SCIP_OKAY;
}

/** output method of statistics table to output file stream 'file' */
static
SCIP_DECL_TABLEOUTPUT(tableOutputPricer)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(table != NULL);

   SCIPprintPricerStatistics(scip, file);

   return SCIP_OKAY;
}

/** output method of statistics table to output file stream 'file' */
static
SCIP_DECL_TABLEOUTPUT(tableOutputBranch)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(table != NULL);

   SCIPprintBranchruleStatistics(scip, file);

   return SCIP_OKAY;
}

/** output method of statistics table to output file stream 'file' */
static
SCIP_DECL_TABLEOUTPUT(tableOutputHeur)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(table != NULL);

   SCIPprintHeuristicStatistics(scip, file);

   return SCIP_OKAY;
}

/** output method of statistics table to output file stream 'file' */
static
SCIP_DECL_TABLEOUTPUT(tableOutputCompression)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(table != NULL);

   SCIPprintCompressionStatistics(scip, file);

   return SCIP_OKAY;
}

/** output method of statistics table to output file stream 'file' */
static
SCIP_DECL_TABLEOUTPUT(tableOutputLP)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(table != NULL);

   SCIPprintLPStatistics(scip, file);

   return SCIP_OKAY;
}

/** output method of statistics table to output file stream 'file' */
static
SCIP_DECL_TABLEOUTPUT(tableOutputNLP)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(table != NULL);

   SCIPprintNLPStatistics(scip, file);

   return SCIP_OKAY;
}

/** output method of statistics table to output file stream 'file' */
static
SCIP_DECL_TABLEOUTPUT(tableOutputRelax)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(table != NULL);

   SCIPprintRelaxatorStatistics(scip, file);

   return SCIP_OKAY;
}

/** output method of statistics table to output file stream 'file' */
static
SCIP_DECL_TABLEOUTPUT(tableOutputTree)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(table != NULL);

   SCIPprintTreeStatistics(scip, file);

   return SCIP_OKAY;
}

/** output method of statistics table to output file stream 'file' */
static
SCIP_DECL_TABLEOUTPUT(tableOutputRoot)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(table != NULL);

   SCIPprintRootStatistics(scip, file);

   return SCIP_OKAY;
}

/** output method of statistics table to output file stream 'file' */
static
SCIP_DECL_TABLEOUTPUT(tableOutputSol)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(table != NULL);

   SCIPprintSolutionStatistics(scip, file);

   return SCIP_OKAY;
}

/** output method of statistics table to output file stream 'file' */
static
SCIP_DECL_TABLEOUTPUT(tableOutputConc)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(table != NULL);

   SCIPprintConcsolverStatistics(scip, file);

   return SCIP_OKAY;
}


/*
 * statistics table specific interface methods
 */

/** creates the default statistics tables and includes them in SCIP */
SCIP_RETCODE SCIPincludeTableDefault(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_TABLE* tmptable;

   tmptable = SCIPfindTable(scip, TABLE_NAME_STATUS);
   if( tmptable == NULL )
   {
      SCIP_CALL( SCIPincludeTable(scip, TABLE_NAME_STATUS, TABLE_DESC_STATUS, SCIP_TABLESTATUS_ON,
            tableCopyDefault, NULL, NULL, NULL, NULL, NULL, tableOutputStatus,
            NULL, TABLE_POSITION_STATUS, TABLE_EARLIEST_STAGE_STATUS) );
   }

   tmptable = SCIPfindTable(scip, TABLE_NAME_TIMING);
   if( tmptable == NULL )
   {
      SCIP_CALL( SCIPincludeTable(scip, TABLE_NAME_TIMING, TABLE_DESC_TIMING, SCIP_TABLESTATUS_ON,
            tableCopyDefault, NULL, NULL, NULL, NULL, NULL, tableOutputTiming,
            NULL, TABLE_POSITION_TIMING, TABLE_EARLIEST_STAGE_TIMING) );
   }

   tmptable = SCIPfindTable(scip, TABLE_NAME_ORIGPROB);
   if( tmptable == NULL )
   {
      SCIP_CALL( SCIPincludeTable(scip, TABLE_NAME_ORIGPROB, TABLE_DESC_ORIGPROB, SCIP_TABLESTATUS_ON,
            tableCopyDefault, NULL, NULL, NULL, NULL, NULL, tableOutputOrigProb,
            NULL, TABLE_POSITION_ORIGPROB, TABLE_EARLIEST_STAGE_ORIGPROB) );
   }

   tmptable = SCIPfindTable(scip, TABLE_NAME_TRANSPROB);
   if( tmptable == NULL )
   {
      SCIP_CALL( SCIPincludeTable(scip, TABLE_NAME_TRANSPROB, TABLE_DESC_TRANSPROB, SCIP_TABLESTATUS_ON,
            tableCopyDefault, NULL, NULL, NULL, NULL, NULL, tableOutputTransProb,
            NULL, TABLE_POSITION_TRANSPROB, TABLE_EARLIEST_STAGE_TRANSPROB) );
   }

   tmptable = SCIPfindTable(scip, TABLE_NAME_PRESOL);
   if( tmptable == NULL )
   {
      SCIP_CALL( SCIPincludeTable(scip, TABLE_NAME_PRESOL, TABLE_DESC_PRESOL, SCIP_TABLESTATUS_ON,
            tableCopyDefault, NULL, NULL, NULL, NULL, NULL, tableOutputPresol,
            NULL, TABLE_POSITION_PRESOL, TABLE_EARLIEST_STAGE_PRESOL) );
   }

   tmptable = SCIPfindTable(scip, TABLE_NAME_CONS);
   if( tmptable == NULL )
   {
      SCIP_CALL( SCIPincludeTable(scip, TABLE_NAME_CONS, TABLE_DESC_CONS, SCIP_TABLESTATUS_ON,
            tableCopyDefault, NULL, NULL, NULL, NULL, NULL, tableOutputCons,
            NULL, TABLE_POSITION_CONS, TABLE_EARLIEST_STAGE_CONS) );
   }

   tmptable = SCIPfindTable(scip, TABLE_NAME_CONSTIMING);
   if( tmptable == NULL )
   {
      SCIP_CALL( SCIPincludeTable(scip, TABLE_NAME_CONSTIMING, TABLE_DESC_CONSTIMING, SCIP_TABLESTATUS_ON,
            tableCopyDefault, NULL, NULL, NULL, NULL, NULL, tableOutputConstiming,
            NULL, TABLE_POSITION_CONSTIMING, TABLE_EARLIEST_STAGE_CONSTIMING) );
   }

   tmptable = SCIPfindTable(scip, TABLE_NAME_PROP);
   if( tmptable == NULL )
   {
      SCIP_CALL( SCIPincludeTable(scip, TABLE_NAME_PROP, TABLE_DESC_PROP, SCIP_TABLESTATUS_ON,
            tableCopyDefault, NULL, NULL, NULL, NULL, NULL, tableOutputProp,
            NULL, TABLE_POSITION_PROP, TABLE_EARLIEST_STAGE_PROP) );
   }

   tmptable = SCIPfindTable(scip, TABLE_NAME_CONFLICT);
   if( tmptable == NULL )
   {
      SCIP_CALL( SCIPincludeTable(scip, TABLE_NAME_CONFLICT, TABLE_DESC_CONFLICT, SCIP_TABLESTATUS_ON,
            tableCopyDefault, NULL, NULL, NULL, NULL, NULL, tableOutputConflict,
            NULL, TABLE_POSITION_CONFLICT, TABLE_EARLIEST_STAGE_CONFLICT) );
   }

   tmptable = SCIPfindTable(scip, TABLE_NAME_SEPA);
   if( tmptable == NULL )
   {
      SCIP_CALL( SCIPincludeTable(scip, TABLE_NAME_SEPA, TABLE_DESC_SEPA, SCIP_TABLESTATUS_ON,
            tableCopyDefault, NULL, NULL, NULL, NULL, NULL, tableOutputSepa,
            NULL, TABLE_POSITION_SEPA, TABLE_EARLIEST_STAGE_SEPA) );
   }

   tmptable = SCIPfindTable(scip, TABLE_NAME_PRICER);
   if( tmptable == NULL )
   {
      SCIP_CALL( SCIPincludeTable(scip, TABLE_NAME_PRICER, TABLE_DESC_PRICER, SCIP_TABLESTATUS_ON,
            tableCopyDefault, NULL, NULL, NULL, NULL, NULL, tableOutputPricer,
            NULL, TABLE_POSITION_PRICER, TABLE_EARLIEST_STAGE_PRICER) );
   }

   tmptable = SCIPfindTable(scip, TABLE_NAME_BRANCH);
   if( tmptable == NULL )
   {
      SCIP_CALL( SCIPincludeTable(scip, TABLE_NAME_BRANCH, TABLE_DESC_BRANCH, SCIP_TABLESTATUS_ON,
            tableCopyDefault, NULL, NULL, NULL, NULL, NULL, tableOutputBranch,
            NULL, TABLE_POSITION_BRANCH, TABLE_EARLIEST_STAGE_BRANCH) );
   }

   tmptable = SCIPfindTable(scip, TABLE_NAME_HEUR);
   if( tmptable == NULL )
   {
      SCIP_CALL( SCIPincludeTable(scip, TABLE_NAME_HEUR, TABLE_DESC_HEUR, SCIP_TABLESTATUS_ON,
            tableCopyDefault, NULL, NULL, NULL, NULL, NULL, tableOutputHeur,
            NULL, TABLE_POSITION_HEUR, TABLE_EARLIEST_STAGE_HEUR) );
   }

   tmptable = SCIPfindTable(scip, TABLE_NAME_COMPRESSION);
   if( tmptable == NULL )
   {
      SCIP_CALL( SCIPincludeTable(scip, TABLE_NAME_COMPRESSION, TABLE_DESC_COMPRESSION, SCIP_TABLESTATUS_ON,
            tableCopyDefault, NULL, NULL, NULL, NULL, NULL, tableOutputCompression,
            NULL, TABLE_POSITION_COMPRESSION, TABLE_EARLIEST_STAGE_COMPRESSION) );
   }

   tmptable = SCIPfindTable(scip, TABLE_NAME_LP);
   if( tmptable == NULL )
   {
      SCIP_CALL( SCIPincludeTable(scip, TABLE_NAME_LP, TABLE_DESC_LP, SCIP_TABLESTATUS_ON,
            tableCopyDefault, NULL, NULL, NULL, NULL, NULL, tableOutputLP,
            NULL, TABLE_POSITION_LP, TABLE_EARLIEST_STAGE_LP) );
   }

   tmptable = SCIPfindTable(scip, TABLE_NAME_NLP);
   if( tmptable == NULL )
   {
      SCIP_CALL( SCIPincludeTable(scip, TABLE_NAME_NLP, TABLE_DESC_NLP, SCIP_TABLESTATUS_ON,
            tableCopyDefault, NULL, NULL, NULL, NULL, NULL, tableOutputNLP,
            NULL, TABLE_POSITION_NLP, TABLE_EARLIEST_STAGE_NLP) );
   }

   tmptable = SCIPfindTable(scip, TABLE_NAME_RELAX);
   if( tmptable == NULL )
   {
      SCIP_CALL( SCIPincludeTable(scip, TABLE_NAME_RELAX, TABLE_DESC_RELAX, SCIP_TABLESTATUS_ON,
            tableCopyDefault, NULL, NULL, NULL, NULL, NULL, tableOutputRelax,
            NULL, TABLE_POSITION_RELAX, TABLE_EARLIEST_STAGE_RELAX) );
   }

   tmptable = SCIPfindTable(scip, TABLE_NAME_TREE);
   if( tmptable == NULL )
   {
      SCIP_CALL( SCIPincludeTable(scip, TABLE_NAME_TREE, TABLE_DESC_TREE, SCIP_TABLESTATUS_ON,
            tableCopyDefault, NULL, NULL, NULL, NULL, NULL, tableOutputTree,
            NULL, TABLE_POSITION_TREE, TABLE_EARLIEST_STAGE_TREE) );
   }

   tmptable = SCIPfindTable(scip, TABLE_NAME_ROOT);
   if( tmptable == NULL )
   {
      SCIP_CALL( SCIPincludeTable(scip, TABLE_NAME_ROOT, TABLE_DESC_ROOT, SCIP_TABLESTATUS_ON,
            tableCopyDefault, NULL, NULL, NULL, NULL, NULL, tableOutputRoot,
            NULL, TABLE_POSITION_ROOT, TABLE_EARLIEST_STAGE_ROOT) );
   }

   tmptable = SCIPfindTable(scip, TABLE_NAME_SOL);
   if( tmptable == NULL )
   {
      SCIP_CALL( SCIPincludeTable(scip, TABLE_NAME_SOL, TABLE_DESC_SOL, SCIP_TABLESTATUS_ON,
            tableCopyDefault, NULL, NULL, NULL, NULL, NULL, tableOutputSol,
            NULL, TABLE_POSITION_SOL, TABLE_EARLIEST_STAGE_SOL) );
   }

   tmptable = SCIPfindTable(scip, TABLE_NAME_CONC);
   if( tmptable == NULL )
   {
      SCIP_CALL( SCIPincludeTable(scip, TABLE_NAME_CONC, TABLE_DESC_CONC, SCIP_TABLESTATUS_ON,
            tableCopyDefault, NULL, NULL, NULL, NULL, NULL, tableOutputConc,
            NULL, TABLE_POSITION_CONC, TABLE_EARLIEST_STAGE_CONC) );
   }

   return SCIP_OKAY;
}
