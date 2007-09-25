/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2007 Tobias Achterberg                              */
/*                                                                           */
/*                  2002-2007 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#pragma ident "@(#) $Id: sepa_mcf.c,v 1.4 2007/09/25 14:32:23 bzfpfend Exp $"

#define SCIP_DEBUG
/**@file   sepa_mcf.c
 * @brief  multi-commodity-flow network cut separator
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>

#include "scip/sepa_mcf.h"


#define SEPA_NAME              "mcf"
#define SEPA_DESC              "multi-commodity-flow network cut separator"
#define SEPA_PRIORITY            -10000
#define SEPA_FREQ                     0
#define SEPA_MAXBOUNDDIST           0.0
#define SEPA_DELAY                FALSE /**< should separation method be delayed, if other separators found cuts? */




/*
 * Data structures
 */

struct SCIP_McfNetwork
{
   SCIP_ROW***           nodeflowconss;      /**< nodeflowconss[v][k]: flow conservation constraint for node v and
                                              *   commodity k; NULL if this is the "missing" constraint */
   SCIP_Real**           nodeflowscales;     /**< scaling factors to convert nodeflowconss[v][k] into a +/-1 <= row */
   SCIP_COL***           arcflows;           /**< arcflows[a][k]: flow variable for arc a and commodity k */
   SCIP_ROW***           arccapacityconss;   /**< arccapacity[a][r]: capacity constraint on arc a for resource r;
                                              *   NULL if uncapacitated */
   SCIP_Real**           arccapacityscales;  /**< scaling factors to convert arccapacity[a][r] into a <= row with
                                              *   positive entries for the flow variables */
   int                   nnodes;             /**< number of nodes in the flow */
   int                   narcs;              /**< number of arcs */
   int                   ncommodities;       /**< number of commodities */
   int                   nresources;         /**< number of resources */
};
typedef struct SCIP_McfNetwork SCIP_MCFNETWORK;


/** separator data */
struct SCIP_SepaData
{
   SCIP_MCFNETWORK*       mcfnetwork;        /**< multi-commodity-flow network structure */
};




/*
 * Local methods
 */

/** creates an empty MCF network data structure */
static
SCIP_RETCODE mcfnetworkCreate(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_MCFNETWORK**     mcfnetwork          /**< MCF network structure */
   )
{
   assert(mcfnetwork != NULL);
   assert(*mcfnetwork == NULL);

   SCIP_CALL( SCIPallocMemory(scip, mcfnetwork) );
   (*mcfnetwork)->nodeflowconss = NULL;
   (*mcfnetwork)->nodeflowscales = NULL;
   (*mcfnetwork)->arcflows = NULL;
   (*mcfnetwork)->arccapacityconss = NULL;
   (*mcfnetwork)->arccapacityscales = NULL;
   (*mcfnetwork)->nnodes = 0;
   (*mcfnetwork)->narcs = 0;
   (*mcfnetwork)->ncommodities = 0;
   (*mcfnetwork)->nresources = 0;

   return SCIP_OKAY;
}

/** frees MCF network data structure */
static
void mcfnetworkFree(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_MCFNETWORK**     mcfnetwork          /**< MCF network structure */
   )
{
   assert(mcfnetwork != NULL);

   if( *mcfnetwork != NULL )
   {
      int v;
      int a;

      for( v = 0; v < (*mcfnetwork)->nnodes; v++ )
      {
         SCIPfreeMemoryArrayNull(scip, &(*mcfnetwork)->nodeflowconss[v]);
         SCIPfreeMemoryArrayNull(scip, &(*mcfnetwork)->nodeflowscales[v]);
      }
      SCIPfreeMemoryArrayNull(scip, &(*mcfnetwork)->nodeflowconss);
      SCIPfreeMemoryArrayNull(scip, &(*mcfnetwork)->nodeflowscales);
      
      for( a = 0; a < (*mcfnetwork)->narcs; a++ )
      {
         SCIPfreeMemoryArrayNull(scip, &(*mcfnetwork)->arcflows[a]);
         SCIPfreeMemoryArrayNull(scip, &(*mcfnetwork)->arccapacityconss[a]);
         SCIPfreeMemoryArrayNull(scip, &(*mcfnetwork)->arccapacityscales[a]);
      }
      SCIPfreeMemoryArrayNull(scip, &(*mcfnetwork)->arcflows);
      SCIPfreeMemoryArrayNull(scip, &(*mcfnetwork)->arccapacityconss);
      SCIPfreeMemoryArrayNull(scip, &(*mcfnetwork)->arccapacityscales);

      SCIPfreeMemory(scip, mcfnetwork);
   }
}


#define LHSPOSSIBLE 1
#define RHSPOSSIBLE 2

/** extracts a MCF network structure from the current LP */
static
SCIP_RETCODE mcfnetworkExtract(
   SCIP*                 scip,               /**< SCIP data structure */ 
   SCIP_MCFNETWORK**     mcfnetwork          /**< MCF network structure */
   )
{
   SCIP_MCFNETWORK* mcf;
   SCIP_ROW** rows;
   SCIP_COL** cols;
   unsigned char* flowrowsigns;
   SCIP_Real* flowrowscalars;
   unsigned char* capacityrowsigns;
   int* pluscom;
   int* minuscom;
   int nrows;
   int ncols;
   int r;

   assert(mcfnetwork != NULL);

   /* create network data structure */
   SCIP_CALL( mcfnetworkCreate(scip, mcfnetwork) );
   assert(*mcfnetwork != NULL);
   mcf = *mcfnetwork;

   /* get LP data */
   SCIP_CALL( SCIPgetLPRowsData(scip, &rows, &nrows) );
   SCIP_CALL( SCIPgetLPColsData(scip, &cols, &ncols) );

   /* extract candidates for flow conservation and capacity constraints */
   SCIP_CALL( SCIPallocBufferArray(scip, &flowrowsigns, nrows) );
   SCIP_CALL( SCIPallocBufferArray(scip, &flowrowscalars, nrows) );
   SCIP_CALL( SCIPallocBufferArray(scip, &capacityrowsigns, nrows) );
   for( r = 0; r < nrows; r++ )
   {
      SCIP_ROW* row;
      SCIP_COL** rowcols;
      SCIP_Real* rowvals;
      SCIP_Real rowlhs;
      SCIP_Real rowrhs;
      int rowlen;
      SCIP_Real coef;
      SCIP_Bool hasinteger;
      unsigned int rowsign;
      int i;

      row = rows[r];
      assert(SCIProwGetLPPos(row) == r);

      flowrowsigns[r] = 0;
      flowrowscalars[r] = 0.0;
      capacityrowsigns[r] = 0;

      rowlen = SCIProwGetNNonz(row);
      if( rowlen == 0 )
         continue;
      rowcols = SCIProwGetCols(row);
      rowvals = SCIProwGetVals(row);
      rowlhs = SCIProwGetLhs(row);
      rowrhs = SCIProwGetRhs(row);

      /* identify flow conservation constraints */
      coef = ABS(rowvals[0]);
      for( i = 1; i < rowlen; i++ )
      {
         SCIP_Real absval = ABS(rowvals[i]);
         if( !SCIPisEQ(scip, absval, coef) )
            break;
      }
      if( i == rowlen )
      {
         if( !SCIPisInfinity(scip, -rowlhs) )
            flowrowsigns[r] |= LHSPOSSIBLE;
         if( !SCIPisInfinity(scip, rowrhs) )
            flowrowsigns[r] |= RHSPOSSIBLE;
         flowrowscalars[r] = 1.0/coef;
      }

      /* identify capacity constraints */
      hasinteger = FALSE;
      rowsign = 0;
      if( !SCIPisInfinity(scip, -rowlhs) )
         rowsign |= LHSPOSSIBLE;
      if( !SCIPisInfinity(scip, rowrhs) )
         rowsign |= RHSPOSSIBLE;
      for( i = 0; i < rowlen && rowsign > 0; i++ )
      {
         if( SCIPvarGetType(SCIPcolGetVar(rowcols[i])) == SCIP_VARTYPE_CONTINUOUS )
         {
            if( rowvals[i] > 0.0 )
               rowsign &= ~LHSPOSSIBLE;
            else
               rowsign &= ~RHSPOSSIBLE;
         }
         else
         {
            hasinteger = TRUE;
            if( SCIPvarGetType(SCIPcolGetVar(rowcols[i])) != SCIP_VARTYPE_BINARY )
            {
               if( rowvals[i] < 0.0 )
                  rowsign &= ~LHSPOSSIBLE;
               else
                  rowsign &= ~RHSPOSSIBLE;
            }
         }
      }
      if( i == rowlen && hasinteger )
         capacityrowsigns[r] = rowsign;
   }

   SCIPdebugMessage("flow conservation candidates:\n");
   for( r = 0; r < nrows; r++ )
   {
      if( flowrowsigns[r] != 0 )
      {
         //SCIPdebug(SCIPprintRow(scip, rows[r], NULL));
         SCIPdebugMessage("%s\n", SCIProwGetName(rows[r]));
      }
   }
   SCIPdebugMessage("capacity candidates:\n");
   for( r = 0; r < nrows; r++ )
   {
      if( capacityrowsigns[r] != 0 )
      {
         //SCIPdebug(SCIPprintRow(scip, rows[r], NULL));
         SCIPdebugMessage("%s\n", SCIProwGetName(rows[r]));
      }
   }

   /* Algorithm to identify multi-commodity-flow network with capacity constraints
    *
    * 1. Initialize pluscom[c] = -1 and minuscom[c] = -1 for all columns c.
    * 2. For all commodities k: Find next node flow conservation constraint for commodity k.
    *    (a) For all columns already known in commodity k with only one entry:
    *        (i) For all flow conservation rows in which this column appears:
    *            - Check if it fits to this commodity. If yes, take it and goto 3.
    *    (b) If nothing fits, exit.
    * 3. For all newly introduced arc flow variables of commodity 0 due to the new row,
    *    find capacity constraints.
    * 4. For first new capacity constraint (just pick one, it is only used as search
    *    criterion for new flow conservation constraints in other commodities):
    *    (a) For all potential arc flow variables in the capacity constraint:
    *        (i) For all flow conservation rows in which this column appears:
    *            - Check if it fits to one of the known commodities.
    *               * plus/minus pattern with overlap
    *               * capacity constraints
    *              If yes, add it to the commodity.
    *            - Otherwise, check if it gives rise to a new commodity.
    *               * plus/minus pattern without overlap
    *               * capacity constraints
    *              If yes, generate new commodity and add the row.
    * 5. Goto 2.
    */
   SCIP_CALL( SCIPallocBufferArray(scip, &pluscom, &ncols) );
   SCIP_CALL( SCIPallocBufferArray(scip, &minuscom, &ncols) );

   /* 1. Initialize pluscom[c] = -1 and minuscom[c] = -1 for all columns c. */
   for( c = 0; c < ncols; c++ )
   {
      pluscom[c] = -1;
      minuscom[c] = -1;
   }

   while( TRUE )
   {
      /* 2. For all commodities k: Find next node flow conservation constraint for commodity k.
       *    (a) For all columns already known in commodity k with only one entry:
       *        (i) For all flow conservation rows in which this column appears:
       *            - Check if it fits to this commodity. If yes, take it and goto 3.
       *    (b) If nothing fits, exit.
       */
      
   }

   /* free memory */
   SCIPfreeBufferArray(scip, &minuscom);
   SCIPfreeBufferArray(scip, &pluscom);
   SCIPfreeBufferArray(scip, &capacityrowsigns);
   SCIPfreeBufferArray(scip, &flowrowscalars);
   SCIPfreeBufferArray(scip, &flowrowsigns);

   return SCIP_OKAY;
}




/*
 * Callback methods of separator
 */

/* TODO: Implement all necessary separator methods. The methods with an #if 0 ... #else #define ... are optional */

/** destructor of separator to free user data (called when SCIP is exiting) */
static
SCIP_DECL_SEPAFREE(sepaFreeMcf)
{  /*lint --e{715}*/
   SCIP_SEPADATA* sepadata;

   /* free separator data */
   sepadata = SCIPsepaGetData(sepa);
   assert(sepadata != NULL);

   /* free MCF network */
   mcfnetworkFree(scip, &sepadata->mcfnetwork);

   SCIPfreeMemory(scip, &sepadata);

   SCIPsepaSetData(sepa, NULL);

   return SCIP_OKAY;
}


/** initialization method of separator (called after problem was transformed) */
#if 0
static
SCIP_DECL_SEPAINIT(sepaInitMcf)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of mcf separator not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define sepaInitMcf NULL
#endif


/** deinitialization method of separator (called before transformed problem is freed) */
#if 0
static
SCIP_DECL_SEPAEXIT(sepaExitMcf)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of mcf separator not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define sepaExitMcf NULL
#endif


/** solving process initialization method of separator (called when branch and bound process is about to begin) */
#if 0
static
SCIP_DECL_SEPAINITSOL(sepaInitsolMcf)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of mcf separator not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define sepaInitsolMcf NULL
#endif


/** solving process deinitialization method of separator (called before branch and bound process data is freed) */
#if 0
static
SCIP_DECL_SEPAEXITSOL(sepaExitsolMcf)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of mcf separator not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define sepaExitsolMcf NULL
#endif


/** LP solution separation method of separator */
static
SCIP_DECL_SEPAEXECLP(sepaExeclpMcf)
{  /*lint --e{715}*/
   SCIP_SEPADATA* sepadata;

   /* get separator data */
   sepadata = SCIPsepaGetData(sepa);
   assert(sepadata != NULL);

   /* get or extract network flow structure */
   if( sepadata->mcfnetwork == NULL )
   {
      SCIP_CALL( mcfnetworkExtract(scip, &sepadata->mcfnetwork) );
      SCIPdebugMessage("extracted network has %d nodes, %d arcs, %d commodities, and %d resources\n",
                       sepadata->mcfnetwork->nnodes, sepadata->mcfnetwork->narcs,
                       sepadata->mcfnetwork->ncommodities, sepadata->mcfnetwork->nresources);
   }
   assert(sepadata->mcfnetwork != NULL);

   SCIPerrorMessage("method of mcf separator not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}


/** arbitrary primal solution separation method of separator */
#if 0
static
SCIP_DECL_SEPAEXECSOL(sepaExecsolMcf)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of mcf separator not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define sepaExecsolMcf NULL
#endif




/*
 * separator specific interface methods
 */

/** creates the mcf separator and includes it in SCIP */
SCIP_RETCODE SCIPincludeSepaMcf(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_SEPADATA* sepadata;

   /* disabled, because separator is not yet finished */
   return SCIP_OKAY;

   /* create cmir separator data */
   SCIP_CALL( SCIPallocMemory(scip, &sepadata) );
   sepadata->mcfnetwork = NULL;

   /* include separator */
   SCIP_CALL( SCIPincludeSepa(scip, SEPA_NAME, SEPA_DESC, SEPA_PRIORITY, SEPA_FREQ, SEPA_MAXBOUNDDIST, SEPA_DELAY,
         sepaFreeMcf, sepaInitMcf, sepaExitMcf, 
         sepaInitsolMcf, sepaExitsolMcf,
         sepaExeclpMcf, sepaExecsolMcf,
         sepadata) );

   /* add mcf separator parameters */
   /* TODO: (optional) add separator specific parameters with SCIPaddTypeParam() here */

   return SCIP_OKAY;
}
