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
#pragma ident "@(#) $Id: sepa_mcf.c,v 1.7 2007/12/10 15:08:17 bzfpfend Exp $"

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

#define LHSPOSSIBLE 1
#define RHSPOSSIBLE 2
#define LHSASSIGNED 4
#define RHSASSIGNED 8

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

/* generates a new commodity */
static
SCIP_RETCODE newCommodity(
   SCIP*                 scip,               /**< SCIP data structure */
   int***                comcols,            /**< pointer to comcols[][] array */
   int**                 ncomcols,           /**< pointer to ncomcols[] array */
   int**                 comcolssize,        /**< pointer to comcolssize[] array */
   int*                  ncommodities,       /**< pointer to number of commodities */
   int*                  commoditiessize     /**< pointer to commoditiessize */
   )
{
   assert(*ncommodities <= *commoditiessize);
   if( *ncommodities == *commoditiessize )
   {
      *commoditiessize *= 2;
      *commoditiessize = MAX(*commoditiessize, (*ncommodities)+1);
      SCIP_CALL( SCIPreallocMemoryArray(scip, comcols, *commoditiessize) );
      SCIP_CALL( SCIPreallocMemoryArray(scip, ncomcols, *commoditiessize) );
      SCIP_CALL( SCIPreallocMemoryArray(scip, comcolssize, *commoditiessize) );
   }
   assert(*ncommodities < *commoditiessize);
   (*ncomcols)[*ncommodities] = 0;
   (*comcolssize)[*ncommodities] = 8;
   SCIP_CALL( SCIPallocMemoryArray(scip, &(*comcols)[*ncommodities], (*comcolssize)[*ncommodities]) );
   (*ncommodities)++;
   
   return SCIP_OKAY;
}

/* add column c to commodity k */
static
SCIP_RETCODE addCommodityCol(
   SCIP*                 scip,               /**< SCIP data structure */
   int**                 comcols,            /**< array to store columns for each commodity */
   int*                  ncomcols,           /**< number of columns in each commodity */
   int*                  comcolssize,        /**< size of comcols[k] arrays */
   int                   k,                  /**< commodity to add column to */
   int                   c                   /**< column to add */
   )
{
   assert(ncomcols[k] <= comcolssize[k]);
   if( ncomcols[k] == comcolssize[k] )
   {
      comcolssize[k] *= 2;
      comcolssize[k] = MAX(comcolssize[k], ncomcols[k]+1);
      SCIP_CALL( SCIPreallocMemoryArray(scip, &comcols[k], comcolssize[k]) );
   }
   assert(ncomcols[k] < comcolssize[k]);
   comcols[k][ncomcols[k]] = c;
   ncomcols[k]++;

   return SCIP_OKAY;
}

/* adds a flow row to the given commodity */
static
SCIP_RETCODE addCommodityFlowRow(
   SCIP*                 scip,               /**< SCIP data structure */
   unsigned char*        flowrowsigns,       /**< marker whether rows are flow row candidates */
   int*                  flowrowcom,         /**< commodities to which rows are assigned */
   int*                  pluscom,            /**< commodities to which columns are assigned with positive entry */
   int*                  minuscom,           /**< commodities to which columns are assigned with negative entry */
   int**                 comcols,            /**< array to store columns for each commodity */
   int*                  ncomcols,           /**< number of columns in each commodity */
   int*                  comcolssize,        /**< size of comcols[k] arrays */
   SCIP_ROW*             row,                /**< flow row to add to commodity */
   unsigned char         flowrowsign,        /**< possible signs to use for the flow row */
   int                   k                   /**< commodity number */
   )
{
   SCIP_COL** rowcols;
   SCIP_Real* rowvals;
   int rowlen;
   int scale;
   int r;
   int i;

   r = SCIProwGetLPPos(row);
   assert(r >= 0);
   assert((flowrowsign & flowrowsigns[r]) == flowrowsign);
   assert((flowrowsign & (LHSPOSSIBLE | RHSPOSSIBLE)) != 0);

   /* prefer lhs
    * TODO: check if we can do something more clever.
    */
   if( (flowrowsign & LHSPOSSIBLE) != 0 )
   {
      flowrowsigns[r] |= LHSASSIGNED;
      scale = +1;
   }
   else
   {
      assert(flowrowsign == RHSPOSSIBLE);
      flowrowsigns[r] |= RHSASSIGNED;
      scale = -1;
   }
   rowcols = SCIProwGetCols(row);
   rowvals = SCIProwGetVals(row);
   rowlen = SCIProwGetNNonz(row);
   for( i = 0; i < rowlen; i++ )
   {
      int c;

      c = SCIPcolGetLPPos(rowcols[i]);
      if( c >= 0 )
      {
         /* if not yet existing, add column to commodity array */
         if( pluscom[c] == -1 && minuscom[c] == -1 )
         {
            SCIP_CALL( addCommodityCol(scip, comcols, ncomcols, comcolssize, k, c) );
         }

         /* assign the sign of the column in the commodity */
         if( scale * rowvals[i] > 0 )
         {
            assert(pluscom[c] == -1);
            pluscom[c] = k;
         }
         else
         {
            assert(minuscom[c] == -1);
            minuscom[c] = k;
         }
      }
   }
   flowrowcom[r] = k;

   return SCIP_OKAY;
}

/* checks whether the given flow constraint candidate row fits into the given commodity, and returns
 * the possible row signs
 */
static
unsigned char getFlowRowCommodityFit(
   unsigned char*        flowrowsigns,       /**< marker whether rows are flow row candidates */
   int*                  flowrowcom,         /**< commodities to which rows are assigned */
   SCIP_Bool*            iscapacity,         /**< marks whether column is identified as capacity variable */
   int*                  pluscom,            /**< commodities to which columns are assigned with positive entry */
   int*                  minuscom,           /**< commodities to which columns are assigned with negative entry */
   SCIP_ROW*             row,                /**< row to check */
   int                   k                   /**< commodity number */
   )
{
   SCIP_COL** rowcols;
   SCIP_Real* rowvals;
   int rowlen;
   unsigned char flowrowsign;
   int r;
   int i;

   r = SCIProwGetLPPos(row);

   /* ignore rows that are not in the current LP */
   if( r < 0 )
      return 0;

   /* check if row is a flow conservation constraint candidate */
   flowrowsign = flowrowsigns[r];
   if( flowrowsign == 0 )
      return 0;

   /* check if row has already been assigned to some commodity */
   if( flowrowcom[r] >= 0 )
      return 0;

   /* flowrowsign must be a combination of the POSSIBLE flags */
   assert((flowrowsign & ~(LHSPOSSIBLE | RHSPOSSIBLE)) == 0);

   /* go through the coefficients and check whether the columns fit into the commodity */
   rowcols = SCIProwGetCols(row);
   rowvals = SCIProwGetVals(row);
   rowlen = SCIProwGetNNonz(row);
   for( i = 0; i < rowlen && flowrowsign != 0; i++ )
   {
      int c;

      c = SCIPcolGetLPPos(rowcols[i]);

      /* ignore columns that are not in the current LP */
      if( c < 0 )
         continue;

      /* capacity variables must not be member of a flow conservation constraint */
      if( iscapacity[c] )
      {
         /* delete row from candidate set */
         flowrowsigns[r] = 0;
         return 0;
      }

      /* if column already belongs to a different commodity, the row cannot be a flow conservation constraint */
      if( (pluscom[c] != k && pluscom[c] != -1)
          || (minuscom[c] != k && minuscom[c] != -1) )
      {
         /* delete row from candidate set */
         flowrowsigns[r] = 0;
         return 0;
      }

      /* if column already belongs to the current commodity, it only fits if it is not yet present with the same sign */
      if( pluscom[c] == k )
      {
         /* column must be included with negative sign */
         if( rowvals[i] < 0.0 )
            flowrowsign &= ~RHSPOSSIBLE;
         else
            flowrowsign &= ~LHSPOSSIBLE;
      }
      if( minuscom[c] == k )
      {
         /* column must be included with positive sign */
         if( rowvals[i] > 0.0 )
            flowrowsign &= ~RHSPOSSIBLE;
         else
            flowrowsign &= ~LHSPOSSIBLE;
      }
   }

   return flowrowsign;
}

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
   int* flowcands;     /* possibly not needed??????????????????? */
   int* capacitycands;
   int nrows;
   int ncols;
   int nflowcands;
   int ncapacitycands;
   int r;
   int c;

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
   SCIP_CALL( SCIPallocBufferArray(scip, &flowcands, nrows) );
   SCIP_CALL( SCIPallocBufferArray(scip, &capacitycands, nrows) );
   nflowcands = 0;
   ncapacitycands = 0;
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
         /* Flow conservation constraints should always be a*x >= d with d >= 0.
          * Depending on the sign of lhs and rhs, we can rule out some options.
          * Only if lhs = 0 and rhs = 0, both sides are still valid candidates.
          */
         if( !SCIPisNegative(scip, rowlhs) )
            flowrowsigns[r] |= LHSPOSSIBLE;
         if( !SCIPisPositive(scip, rowrhs) )
            flowrowsigns[r] |= RHSPOSSIBLE;
         flowrowscalars[r] = 1.0/coef;
         flowcands[nflowcands] = r;
         nflowcands++;
      }

      /* identify capacity constraints */
      hasinteger = FALSE;
      rowsign = 0;
      if( !SCIPisInfinity(scip, -rowlhs) )
         rowsign |= LHSPOSSIBLE;
      if( !SCIPisInfinity(scip, rowrhs) )
         rowsign |= RHSPOSSIBLE;
      for( i = 0; i < rowlen && rowsign != 0; i++ )
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
      if( i == rowlen && hasinteger && rowsign != 0 )
      {
         capacityrowsigns[r] = rowsign;
         capacitycands[ncapacitycands] = r;
         ncapacitycands++;
      }
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

   if( nflowcands > 0 && ncapacitycands > 0 )
   {
      SCIP_Bool* iscapacity; /* marks if column is identified to be a capacity variable */
      int* pluscom;     /* pluscom[c] = k: arc flow column c belongs to commodity k with +1 entry in flow constraint */
      int* minuscom;    /* minuscom[c] = k: arc flow column c belongs to commodity k with -1 entry in flow constraint */
      int* flowrowcom;  /* flowrowcom[r] = k: row r is a flow constraint of commodity k */
      int** comcols;    /* comcols[k]: array of arc flow columns that belong to commodity k */
      int* ncomcols;    /* ncolmcols[k]: number of arc flow columns in commodity k */
      int* comcolssize; /* comcolssize[k]: space available for comcols[] array */
      int ncommodities;
      int commoditiessize;
      int k;

      /* Algorithm to identify multi-commodity-flow network with capacity constraints
       *
       * 1. Sort flow conservation and capacity constraint candidates by a ranking on
       *    how sure we are that it is indeed a constraint of the desired type.
       * 2. Initialize iscapacity[c] = FALSE, pluscom[c] = -1, and minuscom[c] = -1 for all columns c.
       *    Initialize flowrowcom[r] = -1 for all rows r.
       * 3. Put the first flow conservation constraint candidate to commodity 0, and
       *    assign the included variables to commodity 0.
       * 4. Find next node flow conservation constraint (usually in commodity 0).
       *    For all commodities k:
       *     (a) For all columns already known in commodity k with only one entry:
       *          (i) For all flow conservation rows in which this column appears:
       *              - Check if it fits to this commodity. If yes, take it and goto 5.
       *    If nothing fits, exit.
       * 5. For all newly introduced arc flow variables of commodity k due to the new row,
       *    find a capacity constraint which is compatible to all commodities:
       *     (i)   must have entry for arc flow variable
       *     (ii)  must not have entry for other arc flow variables in the commodity k
       *     (iii) must have at most (usually exactly) one entry for arc flow variables of each commodity
       *     (iv)  each entry that belongs to an arc flow variable of a commodity must not have the opposite
       *           sign in this commodity compared to commodity k
       *    If no capacity constraint has been found, goto 4.
       * 6. For the capacity constraint, identify corresponding columns in other commodities:
       *    (a) For all potential arc flow variables in the capacity constraint that do not yet belong to a
       *        commodity:
       *        (i)  For all flow conservation rows in which this column appears:
       *             - Check if it fits to the new node flow conservation constraint in commodity k:
       *                * columns already assigned to a commodity must have the same sign as in new node flow conservation constraint
       *                * number of positive and negative coefficients of unassigned columns must match new node flow conservation constraint
       *               Otherwise, ignore constraint.
       *               Note that if the constraint has overlap with some commodity, it will automatically fit due to the design of the algorithm.
       *             - If found flow conservation constraint has columns which are already assigned to some commodity,
       *               add it to the same commodity and exit this loop.
       *               Otherwise, remember row.
       *        (ii) If no constraint was added to an existing commodity, generate new commodity and add the remembered row.
       * 7. Goto 4.
       */

      /* allocate memory */
      SCIP_CALL( SCIPallocBufferArray(scip, &iscapacity, ncols) );
      SCIP_CALL( SCIPallocBufferArray(scip, &pluscom, ncols) );
      SCIP_CALL( SCIPallocBufferArray(scip, &minuscom, ncols) );
      SCIP_CALL( SCIPallocBufferArray(scip, &flowrowcom, nrows) );

      ncommodities = 0;
      commoditiessize = 8;
      SCIP_CALL( SCIPallocMemoryArray(scip, &comcols, commoditiessize) );
      SCIP_CALL( SCIPallocMemoryArray(scip, &ncomcols, commoditiessize) );
      SCIP_CALL( SCIPallocMemoryArray(scip, &comcolssize, commoditiessize) );

      /* 1. Sort flow conservation and capacity constraint candidates by a ranking on
       *    how sure we are that it is indeed a constraint of the desired type.
       */
      /* TODO: calculate ranking for flow conservation constraint candidates.
       *       probably: give = 0 equations lower priority since they can be used
       *       in both directions and are therefore more flexible. */
      /* TODO: calculate ranking for capacity constraint candidates */

      /* 2. Initialize iscapacity[c] = FALSE, pluscom[c] = -1, and minuscom[c] = -1 for all columns c.
       *    Initialize flowrowcom[r] = -1 for all rows r.
       */
      for( c = 0; c < ncols; c++ )
      {
         iscapacity[c] = FALSE;
         pluscom[c] = -1;
         minuscom[c] = -1;
      }
      for( r = 0; r < nrows; r++ )
         flowrowcom[r] = -1;

      /* 3. Put the first flow conservation constraint candidate to commodity 0, and
       *    assign the included variables to commodity 0.
       */
      SCIP_CALL( newCommodity(scip, &comcols, &ncomcols, &comcolssize, &ncommodities, &commoditiessize) );
      assert(ncommodities == 1);
      r = flowcands[0];
      assert(flowrowsigns[r] != 0);
      SCIP_CALL( addCommodityFlowRow(scip, flowrowsigns, flowrowcom, pluscom, minuscom, comcols, ncomcols, comcolssize,
            rows[r], flowrowsigns[r], 0) );

      /* big extraction loop */
      while( TRUE )
      {
         SCIP_Bool foundfit;

         /* 4. For all commodities k: Find next node flow conservation constraint for commodity k.
          *    (a) For all columns already known in commodity k with only one entry:
          *        (i) For all flow conservation rows in which this column appears:
          *            - Check if it fits to this commodity. If yes, take it and goto 5.
          *    (b) If nothing fits, exit.
          *
          * TODO: We do the extension of the commodities in a greedy fashion w.r.t. the order of the rows
          *       in the column vectors. It may be better to use some predefined priority order.
          */
         foundfit = FALSE;
         for( k = 0; k < ncommodities; k++ )
         {
            int i;

            for( i = 0; i < ncomcols[k]; i++ )
            {
               c = comcols[k][i];
               assert(!iscapacity[c]);
               assert(pluscom[c] == k || minuscom[c] == k);

               /* check if column is missing a +1 or -1 entry in the commodity */
               if( pluscom[c] == -1 || minuscom[c] == -1 )
               {
                  SCIP_ROW** colrows;
                  int collen;
                  int j;

                  /* loop through sparse column vector to find new flow conservation constraints */
                  colrows = SCIPcolGetRows(cols[c]);
                  collen = SCIPcolGetNNonz(cols[c]);
                  for( j = 0; j < collen; j++ )
                  {
                     SCIP_ROW* row;
                     int flowrowsign;

                     row = colrows[j];

                     /* check if row fits into commodity */
                     flowrowsign = getFlowRowCommodityFit(flowrowsigns, flowrowcom, iscapacity, pluscom, minuscom, row, k);
                     if( flowrowsign != 0 )
                     {
                        /* row fits: insert into commodity */
                        SCIP_CALL( addCommodityFlowRow(scip, flowrowsigns, flowrowcom, pluscom, minuscom, comcols, ncomcols, comcolssize,
                              row, flowrowsign, k) );
                        foundfit = TRUE;
                     }
                  }
               }
            }
         }

         /* if we did not find any new flow conservation constraint, exit */
         if( !foundfit )
            break;
      }

#ifdef SCIP_DEBUG
      /* output commodities */
      for( k = 0; k < ncommodities; k++ )
      {
         printf("\n\nCommodity %d:\n\n", k);
         for( r = 0; r < nrows; r++ )
         {
            if( flowrowcom[r] == k )
            {
               int sign;
               if( (flowrowsigns[r] & LHSASSIGNED) != 0 )
                  sign = +1;
               else if ( (flowrowsigns[r] & RHSASSIGNED) != 0 )
                  sign = -1;
               else
                  sign = 0;
               printf("[%+d] ", sign);
               SCIP_CALL( SCIPprintRow(scip, rows[r], NULL) );
            }
         }
      }
#endif

      /* free memory */
      for( k = 0; k < ncommodities; k++ )
         SCIPfreeMemoryArray(scip, &comcols[k]);
      SCIPfreeMemoryArray(scip, &comcols);
      SCIPfreeMemoryArray(scip, &ncomcols);
      SCIPfreeMemoryArray(scip, &comcolssize);

      SCIPfreeBufferArray(scip, &minuscom);
      SCIPfreeBufferArray(scip, &pluscom);
      SCIPfreeBufferArray(scip, &iscapacity);
   }

   /* free memory */
   SCIPfreeBufferArray(scip, &capacitycands);
   SCIPfreeBufferArray(scip, &flowcands);
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

#if 1
   /* disabled, because separator is not yet finished */
   return SCIP_OKAY;
#endif

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
