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
#pragma ident "@(#) $Id: sepa_mcf.c,v 1.8 2008/01/10 16:59:51 bzfpfend Exp $"

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

#define DEFAULT_MAXSIGNDISTANCE       1 /**< maximum Hamming distance of flow conservation constraint sign patterns of the same node */




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
   SCIP_MCFNETWORK*      mcfnetwork;         /**< multi-commodity-flow network structure */
   int                   maxsigndistance;    /**< maximum Hamming distance of flow conservation constraint sign patterns of the same node */
};




/*
 * Local methods
 */

#define LHSPOSSIBLE 1
#define RHSPOSSIBLE 2
#define LHSASSIGNED 4
#define RHSASSIGNED 8

/** applies the commodity merges to the flow rows */
static
void flushCommodityMerges(
   int                   nrows,              /**< number of rows */
   int*                  commerge,           /**< commodity number into which commodity is merged */
   int*                  flowrowcom          /**< commodities to which rows are assigned */
   )
{
   int r;

   for( r = 0; r < nrows; r++ )
   {
      int repk;
      int k;

      repk = flowrowcom[r];
      if( repk < 0 )
         continue;

      /* search representatant */
      while( commerge[repk] != repk )
         repk = commerge[repk];
      
      /* flatten tree */
      k = flowrowcom[r];
      while( commerge[k] != repk )
      {
         int tmpk = commerge[k];
         commerge[k] = repk;
         k = tmpk;
      }

      /* put row into representative commodity */
      flowrowcom[r] = repk;
   }
}

#ifdef SCIP_DEBUG
static
SCIP_RETCODE printCommodities(
   SCIP*                 scip,               /**< SCIP data structure */
   int                   ncommodities,       /**< number of commodities */
   int*                  commerge,           /**< commodity number into which commodity is merged */
   unsigned char*        flowrowsigns,       /**< marker whether rows are flow row candidates */
   int*                  flowrowcom,         /**< commodities to which rows are assigned */
   int*                  colarcid,           /**< arc id for each column, or -1 if not yet assigned */
   int*                  rowarcnodeid,       /**< arc or node id for each capacity or flow conservation constraint, or -1 if not yet assigned */
   int*                  pluscom,            /**< commodities to which columns are assigned with positive entry */
   int*                  minuscom            /**< commodities to which columns are assigned with negative entry */
   )
{
   SCIP_ROW** rows;
   SCIP_COL** cols;
   int nrows;
   int ncols;
   int k;


   SCIP_CALL( SCIPgetLPRowsData(scip, &rows, &nrows) );
   SCIP_CALL( SCIPgetLPColsData(scip, &cols, &ncols) );

   flushCommodityMerges(nrows, commerge, flowrowcom);

   printf("************************************************\n");
   for( k = 0; k < ncommodities; k++ )
   {
      int r;
      int c;

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
            printf("row %2d: [node:%2d] [%+d] ", r, rowarcnodeid[r], sign);
#if 1
            printf("<%s>\n", SCIProwGetName(rows[r]));
#else
            SCIP_CALL( SCIPprintRow(scip, rows[r], NULL) );
#endif
         }
      }

      for( c = 0; c < ncols; c++ )
      {
         if( pluscom[c] == k || minuscom[c] == k )
            printf("col %2d: [arc:%2d] <%s>\n", c, colarcid[c], SCIPvarGetName(SCIPcolGetVar(cols[c])));
      }
   }
   printf("************************************************\n");

   return SCIP_OKAY;
}
#endif

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
   int**                 commerge,           /**< pointer to commerge[] array */
   int*                  ncommodities,       /**< pointer to number of commodities */
   int*                  commoditiessize     /**< pointer to commoditiessize */
   )
{
   assert(*ncommodities <= *commoditiessize);

   SCIPdebugMessage("creating new commodity %d\n", *ncommodities);

   if( *ncommodities == *commoditiessize )
   {
      *commoditiessize *= 2;
      *commoditiessize = MAX(*commoditiessize, (*ncommodities)+1);
      SCIP_CALL( SCIPreallocMemoryArray(scip, comcols, *commoditiessize) );
      SCIP_CALL( SCIPreallocMemoryArray(scip, ncomcols, *commoditiessize) );
      SCIP_CALL( SCIPreallocMemoryArray(scip, comcolssize, *commoditiessize) );
      SCIP_CALL( SCIPreallocMemoryArray(scip, commerge, *commoditiessize) );
   }
   assert(*ncommodities < *commoditiessize);
   (*ncomcols)[*ncommodities] = 0;
   (*comcolssize)[*ncommodities] = 8;
   (*commerge)[*ncommodities] = *ncommodities;
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

   SCIPdebugMessage("adding column %d to commodity %d\n", c, k);

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

/* merges commodity k1 into commodity k2 */
static
SCIP_RETCODE mergeCommodities(
   SCIP*                 scip,               /**< SCIP data structure */
   int*                  pluscom,            /**< commodities to which columns are assigned with positive entry */
   int*                  minuscom,           /**< commodities to which columns are assigned with negative entry */
   int**                 comcols,            /**< array to store columns for each commodity */
   int*                  ncomcols,           /**< number of columns in each commodity */
   int*                  comcolssize,        /**< size of comcols[k] arrays */
   int*                  commerge,           /**< commodity number into which commodity is merged */
   int                   k1,                 /**< commodity to merge into k2 */
   int                   k2                  /**< commodity that is extended */
   )
{
   int i;

   assert(k1 >= 0);
   assert(k2 >= 0);
   assert(commerge[k1] == k1);
   assert(commerge[k2] == k2);
   assert(ncomcols[k1] <= comcolssize[k1]);
   assert(ncomcols[k2] <= comcolssize[k2]);

   SCIPdebugMessage("merging commodity %d into commodity %d\n", k1, k2);

   if( ncomcols[k2] + ncomcols[k1] > comcolssize[k2] )
   {
      comcolssize[k2] *= 2;
      comcolssize[k2] = MAX(comcolssize[k2], ncomcols[k2]+ncomcols[k1]);
      SCIP_CALL( SCIPreallocMemoryArray(scip, &comcols[k2], comcolssize[k2]) );
   }
   assert(ncomcols[k2] + ncomcols[k1] <= comcolssize[k2]);

   for( i = 0; i < ncomcols[k1]; i++ )
   {
      int c;

      c = comcols[k1][i];
      assert(c >= 0);

      /* update column -> commodity sign assignment */
      assert(pluscom[c] == -1 || pluscom[c] == k1);
      assert(minuscom[c] == -1 || minuscom[c] == k1);
      if( pluscom[c] == k1 )
         pluscom[c] = k2;
      if( minuscom[c] == k1 )
         minuscom[c] = k2;

      /* update commodity -> colums array */
      comcols[k2][ncomcols[k2]] = c;
      ncomcols[k2]++;
   }
   ncomcols[k1] = 0;
   commerge[k1] = k2;

   return SCIP_OKAY;
}

/* adds a flow row to the given commodity;
 * merges commodities if necessary
 */
static
SCIP_RETCODE addCommodityFlowRow(
   SCIP*                 scip,               /**< SCIP data structure */
   unsigned char*        flowrowsigns,       /**< marker whether rows are flow row candidates */
   int*                  flowrowcom,         /**< commodities to which rows are assigned */
   int*                  rowarcnodeid,       /**< arc or node id for each capacity or flow conservation constraint, or -1 if not yet assigned */
   unsigned char*        capacityrowsigns,   /**< marker whether rows are capacity row candidates */
   int*                  pluscom,            /**< commodities to which columns are assigned with positive entry */
   int*                  minuscom,           /**< commodities to which columns are assigned with negative entry */
   int**                 comcols,            /**< array to store columns for each commodity */
   int*                  ncomcols,           /**< number of columns in each commodity */
   int*                  comcolssize,        /**< size of comcols[k] arrays */
   int*                  commerge,           /**< commodity number into which commodity is merged */
   SCIP_ROW*             row,                /**< flow row to add to commodity */
   unsigned char         flowrowsign,        /**< possible signs to use for the flow row */
   int                   k,                  /**< commodity number */
   int                   nodeid              /**< node id to which the flow row belongs */
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

   SCIPdebugMessage("adding flow row %d <%s> with sign %+d to commodity %d\n", r, SCIProwGetName(row), scale, k);
   SCIP_CALL( SCIPprintRow(scip, row, NULL) );

   rowcols = SCIProwGetCols(row);
   rowvals = SCIProwGetVals(row);
   rowlen = SCIProwGetNNonz(row);
   for( i = 0; i < rowlen; i++ )
   {
      int c;

      c = SCIPcolGetLPPos(rowcols[i]);
      if( c >= 0 )
      {
         /* if not yet existing, add column to commodity array;
          * if column belongs to a different commodity, merge commodities
          */
         assert(pluscom[c] == -1 || minuscom[c] == -1 || pluscom[c] == minuscom[c]);
         if( pluscom[c] == -1 && minuscom[c] == -1 )
         {
            SCIP_CALL( addCommodityCol(scip, comcols, ncomcols, comcolssize, k, c) );
         }
         else if( pluscom[c] != -1 && pluscom[c] != k )
         {
            SCIP_CALL( mergeCommodities(scip, pluscom, minuscom, comcols, ncomcols, comcolssize, commerge, pluscom[c], k) );
         }
         else if( minuscom[c] != -1 && minuscom[c] != k )
         {
            SCIP_CALL( mergeCommodities(scip, pluscom, minuscom, comcols, ncomcols, comcolssize, commerge, minuscom[c], k) );
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
   capacityrowsigns[r] = 0;
   rowarcnodeid[r] = nodeid;

   return SCIP_OKAY;
}

/* checks whether the given flow constraint candidate row fits into the given commodity, and returns
 * the possible row signs;
 * if commodities need to be merged in order to make the row fit, it stores the number of the last commodity
 * in k
 */
static
unsigned char getFlowRowCommodityFit(
   unsigned char*        flowrowsigns,       /**< marker whether rows are flow row candidates */
   int*                  flowrowcom,         /**< commodities to which rows are assigned */
   int*                  colarcid,           /**< arc id for each column, or -1 if not yet assigned */
   SCIP_Bool*            iscapacity,         /**< marks whether column is identified as capacity variable */
   int*                  pluscom,            /**< commodities to which columns are assigned with positive entry */
   int*                  minuscom,           /**< commodities to which columns are assigned with negative entry */
   int                   nplussigns,         /**< number of plus signs required in flow conservation constraint, or -1 for no restriction */
   int                   nminussigns,        /**< number of minus signs required in flow conservation constraint, or -1 for no restriction */
   int                   nunknownplussigns,  /**< number of plus signs in reference flow constraint without assigned arc */
   int                   nunknownminussigns, /**< number of minus signs in reference flow constraint without assigned arc */
   int*                  arcpattern,         /**< sign pattern for known arc ids required in flow conservation constraint, or NULL */
   int                   fixedc,             /**< column for which the sign is fixed, or -1 */
   int                   fixedcsign,         /**< fixed sign of the column */
   SCIP_ROW*             row,                /**< row to check */
   int*                  k,                  /**< pointer to commodity number; if value is -1 on input, it will be assigned automatically */
   int*                  signdistance        /**< pointer to store the hamming distance to the required sign pattern, or NULL */
   )
{
   SCIP_COL** rowcols;
   SCIP_Real* rowvals;
   int rowlen;
   int nrowposcoefs;
   int nrownegcoefs;
   int nrowunknownposcoefs;
   int nrowunknownnegcoefs;
   unsigned char flowrowsign;
   int r;
   int i;

   assert(k != NULL);
   assert(fixedcsign == -1 || fixedcsign == +1);

   if( signdistance != NULL )
      *signdistance = INT_MAX;

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
   nrowposcoefs = 0;
   nrownegcoefs = 0;
   nrowunknownposcoefs = 0;
   nrowunknownnegcoefs = 0;
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

      /* if this column is the one of fixed sign, check its sign */
      if( c == fixedc )
      {
         if( rowvals[i] > 0.0 )
         {
            if( fixedcsign == +1 )
               flowrowsign &= ~RHSPOSSIBLE;
            else
               flowrowsign &= ~LHSPOSSIBLE;
         }
         else
         {
            if( fixedcsign == -1 )
               flowrowsign &= ~RHSPOSSIBLE;
            else
               flowrowsign &= ~LHSPOSSIBLE;
         }
      }

      /* if column already belongs to a different commodity, the row cannot be a flow conservation constraint */
      assert(pluscom[c] == -1 || minuscom[c] == -1 || pluscom[c] == minuscom[c]);
      *k = MAX(*k, pluscom[c]);
      *k = MAX(*k, minuscom[c]);

      /* column only fits if it is not yet present with the same sign */
      if( pluscom[c] >= 0 )
      {
         /* column must be included with negative sign */
         if( rowvals[i] < 0.0 )
            flowrowsign &= ~RHSPOSSIBLE;
         else
            flowrowsign &= ~LHSPOSSIBLE;
      }
      if( minuscom[c] >= 0 )
      {
         /* column must be included with positive sign */
         if( rowvals[i] > 0.0 )
            flowrowsign &= ~RHSPOSSIBLE;
         else
            flowrowsign &= ~LHSPOSSIBLE;
      }

      /* count positive and negative coefficients */
      if( rowvals[i] > 0.0 )
         nrowposcoefs++;
      else
         nrownegcoefs++;

      /* compare with predefined arc pattern */
      if( arcpattern != NULL )
      {
         int arcid;

         arcid = colarcid[c];
         if( arcid >= 0 )
         {
            SCIPdebugMessage("  -> arcid=%d arcpattern=%+d val=%+g\n", arcid, arcpattern[arcid], rowvals[i]);

            /* check if row has opposite sign for this arc */
            if( arcpattern[arcid] == +1 )
            {
               if( rowvals[i] < 0.0 )
                  flowrowsign &= ~LHSPOSSIBLE;
               else
                  flowrowsign &= ~RHSPOSSIBLE;
            }
            else if( arcpattern[arcid] == -1 )
            {
               if( rowvals[i] > 0.0 )
                  flowrowsign &= ~LHSPOSSIBLE;
               else
                  flowrowsign &= ~RHSPOSSIBLE;
            }
            else
            {
               /* the arc has an id in the current commodity but does not exist in the reference commodity: this cannot be */
               flowrowsign = 0;
            }
         }
         else
         {
            /* flow variable is not assigned to an arc in current commodity: count the number of such variables */
            if( rowvals[i] > 0.0 )
               nrowunknownposcoefs++;
            else
               nrowunknownnegcoefs++;
         }
      }
   }

   /* calculate Hamming distance to required sign pattern */
   if( flowrowsign != 0 && signdistance != NULL )
   {
      int lhsdistance;
      int rhsdistance;

      SCIPdebugMessage("  -> <%s> nsigns: +%d -%d (+%d -%d)  nrowcoefs: +%d -%d (+%d -%d)\n", SCIProwGetName(row), 
                       nplussigns, nminussigns, nunknownplussigns, nunknownminussigns,
                       nrowposcoefs, nrownegcoefs, nrowunknownposcoefs, nrowunknownnegcoefs);
      lhsdistance = 0;
      rhsdistance = 0;
      if( nplussigns >= 0 )
      {
         lhsdistance += ABS(nplussigns - nrowposcoefs);
         rhsdistance += ABS(nplussigns - nrownegcoefs);
      }
      if( nminussigns >= 0 )
      {
         lhsdistance += ABS(nminussigns - nrownegcoefs);
         rhsdistance += ABS(nminussigns - nrowposcoefs);
      }
      if( nunknownplussigns >= 0 )
      {
         lhsdistance += ABS(nunknownplussigns - nrowunknownposcoefs);
         rhsdistance += ABS(nunknownplussigns - nrowunknownnegcoefs);
      }
      if( nunknownminussigns >= 0 )
      {
         lhsdistance += ABS(nunknownminussigns - nrowunknownnegcoefs);
         rhsdistance += ABS(nunknownminussigns - nrowunknownposcoefs);
      }
      if( flowrowsign == LHSPOSSIBLE )
         *signdistance = lhsdistance;
      else if( flowrowsign == RHSPOSSIBLE )
         *signdistance = rhsdistance;
      else if( flowrowsign == (LHSPOSSIBLE | RHSPOSSIBLE) )
      {
         if( lhsdistance < rhsdistance )
            flowrowsign = LHSPOSSIBLE;
         else if( lhsdistance > rhsdistance )
            flowrowsign = RHSPOSSIBLE;
         *signdistance = MIN(lhsdistance, rhsdistance);
      }
   }

   return flowrowsign;
}

/* returns whether the given row can be used as capacity constraint by checking
 *  (i)   must have entry for arc flow variable
 *  (ii)  must not have entry for other arc flow variables in the commodity k
 *  (iii) must have at most (usually exactly) one entry for arc flow variables of each commodity
 *  (iv)  each entry that belongs to an arc flow variable of a commodity must not have the opposite
 *        sign in this commodity compared to commodity k
 */
static
int getCapacityRowFit(
   SCIP*                 scip,               /**< SCIP data structure */ 
   int                   ncommodities,       /**< current number of commodities */
   unsigned char*        capacityrowsigns,   /**< marker whether rows are capacity row candidates */
   int*                  pluscom,            /**< commodities to which columns are assigned with positive entry */
   int*                  minuscom,           /**< commodities to which columns are assigned with negative entry */
   SCIP_Bool             flowsignisplus,     /**< does the new arc flow variable appear with positive sign in new flow conservation constraint? */
   SCIP_ROW*             row                 /**< capacity row candidate */
   )
{
   SCIP_COL** rowcols;
   SCIP_Real* rowvals;
   int rowlen;
   SCIP_Bool* comcovered;
   int* forbidrhssigncom;
   int* forbidlhssigncom;
   int capacityrowsign;
   int r;
   int i;

   r = SCIProwGetLPPos(row);

   /* if row is not in LP, it is no candidate */
   if( r < 0 )
      return 0;

   /* check if row is capacity row candidate */
   capacityrowsign = capacityrowsigns[r];
   if( capacityrowsign == 0 )
      return 0;

   /* initialize commodity coverage flags */
   SCIP_CALL( SCIPallocBufferArray(scip, &comcovered, ncommodities) );
   BMSclearMemoryArray(comcovered, ncommodities);

   /* if entry in new flow conservation constraint is positive, we have to forbid negative signs for arc flow variables */
   if( flowsignisplus )
   {
      forbidrhssigncom = minuscom;
      forbidlhssigncom = pluscom;
   }
   else
   {
      forbidrhssigncom = pluscom;
      forbidlhssigncom = minuscom;
   }

   /* loop through row */
   rowcols = SCIProwGetCols(row);
   rowvals = SCIProwGetVals(row);
   rowlen = SCIProwGetNNonz(row);
   for( i = 0; i < rowlen && capacityrowsign != 0; i++ )
   {
      int c;
      int k;

      c = SCIPcolGetLPPos(rowcols[i]);
      if( c < 0 )
         continue;

      /*  (i)   must have entry for arc flow variable
       *  (ii)  must not have entry for other arc flow variables in the commodity k
       *  (iii) must have at most (usually exactly) one entry for arc flow variables of each commodity
       */
      k = pluscom[c];
      if( k == -1 )
         k = minuscom[c];
      if( k >= 0 )
      {
         if( comcovered[k] )
         {
            capacityrowsign = 0;
            break;
         }
         comcovered[k] = TRUE;
      }

      /*  (iv)  each entry that belongs to an arc flow variable of a commodity must not have the opposite
       *        sign in this commodity compared to commodity k
       */
      if( forbidlhssigncom[c] >= 0 )
         capacityrowsign &= ~LHSPOSSIBLE;
      if( forbidrhssigncom[c] >= 0 )
         capacityrowsign &= ~RHSPOSSIBLE;
   }

   /* free temporary memory */
   SCIPfreeBufferArray(scip, &comcovered);

   return capacityrowsign;
}

/* adds new capacity row to the network and finds flow conservation constraints for variables that are not yet
 * assigned to a commodity:
 *  (i) For all potential arc flow variables in the capacity constraint that do not yet belong to a
 *       commodity:
 *        (1) For all flow conservation rows in which this column appears:
 *            - Check if it fits to the new node flow conservation constraint in commodity k:
 *               * columns already assigned to a commodity must have the same sign as in new node flow conservation constraint
 *               * number of positive and negative coefficients of unassigned columns must match new node flow conservation constraint
 *              Otherwise, ignore constraint.
 *              Note that if the constraint has overlap with some commodity, it will automatically fit due to the design of the algorithm.
 *            - If found flow conservation constraint has columns which are already assigned to some commodity,
 *              add it to the same commodity and exit this loop.
 *              Otherwise, remember row.
 *        (2) If no constraint was added to an existing commodity, generate new commodity and add the remembered row.
 */
static
SCIP_RETCODE addCapacityRow(
   SCIP*                 scip,               /**< SCIP data structure */ 
   int                   maxsigndistance,    /**< maximum Hamming distance of flow conservation constraint sign patterns of the same node */
   unsigned char*        flowrowsigns,       /**< marker whether rows are flow row candidates */
   int*                  flowrowcom,         /**< commodities to which rows are assigned */
   int*                  colarcid,           /**< arc id for each column, or -1 if not yet assigned */
   int*                  rowarcnodeid,       /**< arc or node id for each capacity or flow conservation constraint, or -1 if not yet assigned */
   unsigned char*        capacityrowsigns,   /**< marker whether rows are capacity row candidates */
   SCIP_Bool*            iscapacity,         /**< marks whether column is identified as capacity variable */
   int*                  pluscom,            /**< commodities to which columns are assigned with positive entry */
   int*                  minuscom,           /**< commodities to which columns are assigned with negative entry */
   int                   newrowscale,        /**< scaling (+1 or -1) for new flow conservation constraint */
   int                   nplussigns,         /**< number of plus signs in new flow conservation constraint */
   int                   nminussigns,        /**< number of minus signs in new flow conservation constraint */
   int                   nunknownplussigns,  /**< number of plus signs in reference flow constraint without assigned arc */
   int                   nunknownminussigns, /**< number of minus signs in reference flow constraint without assigned arc */
   int*                  arcpattern,         /**< sign pattern for known arc ids in new flow conservation constraint */
   int***                comcols,            /**< pointer to comcols[][] array */
   int**                 ncomcols,           /**< pointer to ncomcols[] array */
   int**                 comcolssize,        /**< pointer to comcolssize[] array */
   int**                 commerge,           /**< pointer to commerge[] array */
   int*                  ncommodities,       /**< pointer to number of commodities */
   int*                  commoditiessize,    /**< pointer to commoditiessize */
   int                   sourcec,            /**< source column number for which this capacity row was generated */
   SCIP_ROW*             row,                /**< capacity row to add */
   int                   capacityrowsign,    /**< sign to use for the capacity row */
   int                   nodeid,             /**< node id on which we are currently working */
   int                   arcid               /**< arc id to which this capacity row belongs */
   )
{
   SCIP_COL** rowcols;
   int rowlen;
   int sourcecsign;
   int r;
   int i;

   r = SCIProwGetLPPos(row);
   assert(r >= 0);
   assert((flowrowsigns[r] & (LHSASSIGNED | RHSASSIGNED)) == 0);
   assert((capacityrowsigns[r] & capacityrowsign) != 0);
   assert(rowarcnodeid[r] == -1);
   assert(newrowscale == +1 || newrowscale == -1);

   /* sign must have been chosen uniquely */
   if( capacityrowsign == LHSPOSSIBLE )
      capacityrowsigns[r] |= LHSASSIGNED;
   else
   {
      assert(capacityrowsign == RHSPOSSIBLE);
      capacityrowsigns[r] |= RHSASSIGNED;
   }
   flowrowsigns[r] = 0;
   rowarcnodeid[r] = arcid;

   SCIPdebugMessage("adding capacity row %d <%s> with sign %+d\n", r, SCIProwGetName(row),
                    (capacityrowsigns[r] & LHSASSIGNED) != 0 ? -1 : +1);
   SCIP_CALL( SCIPprintRow(scip, row, NULL) );

   /* get the sign of the source column in the flow conservation constraint for which we want to add the capacity constraint;
    * because the column was newly introduced by the flow constraint, its sign can be found in the pluscom/minuscom arrays
    */
   assert((pluscom[sourcec] >= 0) + (minuscom[sourcec] >= 0) == 1);
   if( pluscom[sourcec] >= 0 )
      sourcecsign = +1;
   else
      sourcecsign = -1;

#ifdef SCIP_DEBUG
   {
      int a;
      SCIPdebugMessage("  -> arc pattern:");
      for( a = 0; a <= arcid; a++ )
      {
         if( arcpattern[a] != 0 )
            printf(" %d:%+d", a, arcpattern[a]);
      }
      printf(" [unknown: +%d -%d]\n", nunknownplussigns, nunknownminussigns);
   }   
#endif

   /* find flow conservation constraints for included columns */
   rowcols = SCIProwGetCols(row);
   rowlen = SCIProwGetNNonz(row);
   for( i = 0; i < rowlen; i++ )
   {
      SCIP_ROW** colrows;
      int collen;
      unsigned char bestflowrowsign;
      int bestj;
      int bestk;
      int bestsigndistance;
      int bestscale;
      int j;
      int c;

      c = SCIPcolGetLPPos(rowcols[i]);
      if( c < 0 )
         continue;

      assert(colarcid[c] == -1 || colarcid[c] == arcid);
      colarcid[c] = arcid;
      SCIPdebugMessage(" -> assigning column %d <%s> to arc %d\n", c, SCIPvarGetName(SCIPcolGetVar(rowcols[i])), arcid);

      /* if this is the column for which the capacity row was generated, the corresponding flow row
       * already exists
       */
      if( c == sourcec )
         continue;

      /* check if column already belongs to a commodity */
      if( pluscom[c] >= 0 || minuscom[c] >= 0 )
         continue;

      /* if variable is already declared as capacity variable, ignore it */
      if( iscapacity[c] )
         continue;

      /* find flow conservation constraint */
      SCIPdebugMessage(" -> searching flow row for column %d <%s>\n", c, SCIPvarGetName(SCIPcolGetVar(rowcols[i])));
      colrows = SCIPcolGetRows(rowcols[i]);
      collen = SCIPcolGetNNonz(rowcols[i]);
      bestj = -1;
      bestk = -1;
      bestflowrowsign = 0;
      bestsigndistance = maxsigndistance;
      for( j = 0; j < collen; j++ )
      {
         unsigned char flowrowsign;
         int k;
         int signdistance;

         /* - Check if it fits to the new node flow conservation constraint in commodity k:
          *    * columns already assigned to a commodity must have the same sign as in new node flow conservation constraint
          *    * number of positive and negative coefficients of unassigned columns must match new node flow conservation constraint
          */
         k = -1;
         flowrowsign = getFlowRowCommodityFit(flowrowsigns, flowrowcom, colarcid, iscapacity, pluscom, minuscom,
                                              nplussigns, nminussigns, nunknownplussigns, nunknownminussigns, arcpattern, c, sourcecsign,
                                              colrows[j], &k, &signdistance);
         if( flowrowsign != 0 && signdistance <= bestsigndistance )
         {
            int scale;

            if( flowrowsign == LHSPOSSIBLE )
               scale = +1;
            else if( flowrowsign == RHSPOSSIBLE )
               scale = -1;
            else
            {
               assert(flowrowsign == (LHSPOSSIBLE | RHSPOSSIBLE));
               scale = newrowscale;
               if( scale == +1 )
                  flowrowsign = LHSPOSSIBLE;
               else
               {
                  assert(scale == -1);
                  flowrowsign = RHSPOSSIBLE;
               }
            }

            if( signdistance < bestsigndistance || scale == newrowscale || (bestscale != newrowscale && bestk == -1) )
            {
               bestj = j;
               bestk = k;
               bestflowrowsign = flowrowsign;
               bestsigndistance = signdistance;
               bestscale = scale;
               if( signdistance == 0 && k != -1 )
                  break;
            }
         }
      }

      if( bestj != -1 )
      {
         /* check if we need to create a new commodity */
         if( bestk == -1 )
         {
            /* create new commodity */
            SCIP_CALL( newCommodity(scip, comcols, ncomcols, comcolssize, commerge, ncommodities, commoditiessize) );
            bestk = *ncommodities - 1;
         }

         /* add the flow conservation constraint to commodity bestk */
         SCIPdebugMessage(" -> found flow row %d for column %d <%s> in commodity %d with sign distance %d\n",
                          SCIProwGetLPPos(colrows[bestj]), c, SCIPvarGetName(SCIPcolGetVar(rowcols[i])), bestk, bestsigndistance);
         SCIP_CALL( addCommodityFlowRow(scip, flowrowsigns, flowrowcom, rowarcnodeid, capacityrowsigns, pluscom, minuscom,
                                        *comcols, *ncomcols, *comcolssize, *commerge, colrows[bestj], bestflowrowsign, bestk, nodeid) );
      }
      else
      {
         /* if we did not find a flow conservation constraint for this variable, it must be a capacity variable */
         iscapacity[c] = TRUE;
      }
   }

   return SCIP_OKAY;
}


/** extracts a MCF network structure from the current LP */
static
SCIP_RETCODE mcfnetworkExtract(
   SCIP*                 scip,               /**< SCIP data structure */ 
   int                   maxsigndistance,    /**< maximum Hamming distance of flow conservation constraint sign patterns of the same node */
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
         /* Flow conservation constraints should always be a*x >= d.
          * If lhs and rhs are finite, both sides are still valid candidates.
          */
         if( !SCIPisInfinity(scip, -rowlhs) )
            flowrowsigns[r] |= LHSPOSSIBLE;
         if( !SCIPisInfinity(scip, rowrhs) )
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
         SCIPdebugMessage("%4d: %s\n", r, SCIProwGetName(rows[r]));
      }
   }
   SCIPdebugMessage("capacity candidates:\n");
   for( r = 0; r < nrows; r++ )
   {
      if( capacityrowsigns[r] != 0 )
      {
         //SCIPdebug(SCIPprintRow(scip, rows[r], NULL));
         SCIPdebugMessage("%4d: %s\n", r, SCIProwGetName(rows[r]));
      }
   }

   if( nflowcands > 0 && ncapacitycands > 0 )
   {
      SCIP_Bool* iscapacity; /* marks if column is identified to be a capacity variable */
      int* pluscom;      /* pluscom[c] = k: arc flow column c belongs to commodity k with +1 entry in flow constraint */
      int* minuscom;     /* minuscom[c] = k: arc flow column c belongs to commodity k with -1 entry in flow constraint */
      int* flowrowcom;   /* flowrowcom[r] = k: row r is a flow constraint of commodity k */
      int* colarcid;     /* arc id for each column, or -1 if not yet assigned */
      int* rowarcnodeid; /* arc or node id for each capacity or flow conservation constraint, or -1 if not yet assigned */
      int* arcpattern;   /* temporary array to store the partial arc +/- pattern of a node */
      int** comcols;     /* comcols[k]: array of arc flow columns that belong to commodity k */
      int* ncomcols;     /* ncolmcols[k]: number of arc flow columns in commodity k */
      int* comcolssize;  /* comcolssize[k]: space available for comcols[] array */
      int* commerge;     /* commodity number into which commodity has been merged */
      int ncommodities;
      int commoditiessize;
      SCIP_Bool foundfit;
      SCIP_ROW* newrow;
      int narcs;
      int nnodes;
      int k;

      /* Algorithm to identify multi-commodity-flow network with capacity constraints
       *
       * 1. Sort flow conservation and capacity constraint candidates by a ranking on
       *    how sure we are that it is indeed a constraint of the desired type.
       * 2. Initialize iscapacity[c] = FALSE, pluscom[c] = -1, and minuscom[c] = -1 for all columns c.
       *    Initialize flowrowcom[r] = -1 for all rows r.
       * 3. Put the first flow conservation constraint candidate to commodity 0, and
       *    assign the included variables to commodity 0.
       * 4. For all newly introduced arc flow variables of commodity k due to the new row:
       *    (a) find a capacity constraint which is compatible to all commodities:
       *         (i)   must have entry for arc flow variable
       *         (ii)  must not have entry for other arc flow variables in the commodity k
       *         (iii) must have at most (usually exactly) one entry for arc flow variables of each commodity
       *         (iv)  each entry that belongs to an arc flow variable of a commodity must not have the opposite
       *               sign in this commodity compared to commodity k
       *        If no capacity constraint has been found, goto 5.
       *    (b) For the capacity constraint, identify corresponding columns in other commodities:
       *         (i) For all potential arc flow variables in the capacity constraint that do not yet belong to a
       *             commodity:
       *              (1) For all flow conservation rows in which this column appears:
       *                  - Check if it fits to the new node flow conservation constraint in commodity k:
       *                     * columns already assigned to a commodity must have the same sign as in new node flow conservation constraint
       *                     * number of positive and negative coefficients of unassigned columns must match new node flow conservation constraint
       *                    Otherwise, ignore constraint.
       *                    Note that if the constraint has overlap with some commodity, it will automatically fit due to the design of the algorithm.
       *                  - If found flow conservation constraint has columns which are already assigned to some commodity,
       *                    add it to the same commodity and exit this loop.
       *                    Otherwise, remember row.
       *              (2) If no constraint was added to an existing commodity, generate new commodity and add the remembered row.
       * 5. Find next node flow conservation constraint (usually in commodity 0).
       *    For all commodities k:
       *     (a) For all columns already known in commodity k with only one entry:
       *          (i) For all flow conservation rows in which this column appears:
       *              - Check if it fits to this commodity. If yes, take it and goto 4.
       *    If nothing fits, exit.
       */

      /* allocate memory */
      SCIP_CALL( SCIPallocBufferArray(scip, &iscapacity, ncols) );
      SCIP_CALL( SCIPallocBufferArray(scip, &pluscom, ncols) );
      SCIP_CALL( SCIPallocBufferArray(scip, &minuscom, ncols) );
      SCIP_CALL( SCIPallocBufferArray(scip, &flowrowcom, nrows) );
      SCIP_CALL( SCIPallocBufferArray(scip, &colarcid, ncols) );
      SCIP_CALL( SCIPallocBufferArray(scip, &rowarcnodeid, nrows) );
      SCIP_CALL( SCIPallocBufferArray(scip, &arcpattern, ncols) );

      ncommodities = 0;
      commoditiessize = 8;
      SCIP_CALL( SCIPallocMemoryArray(scip, &comcols, commoditiessize) );
      SCIP_CALL( SCIPallocMemoryArray(scip, &ncomcols, commoditiessize) );
      SCIP_CALL( SCIPallocMemoryArray(scip, &comcolssize, commoditiessize) );
      SCIP_CALL( SCIPallocMemoryArray(scip, &commerge, commoditiessize) );

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
         colarcid[c] = -1;
         arcpattern[c] = 0;
      }
      for( r = 0; r < nrows; r++ )
      {
         flowrowcom[r] = -1;
         rowarcnodeid[r] = -1;
      }
      narcs = 0;
      nnodes = 0;

      SCIPdebugMessage("starting main network extraction\n");

      /* 3. Put the first flow conservation constraint candidate to commodity 0, and
       *    assign the included variables to commodity 0.
       */
      SCIP_CALL( newCommodity(scip, &comcols, &ncomcols, &comcolssize, &commerge, &ncommodities, &commoditiessize) );
      assert(ncommodities == 1);
      r = flowcands[0];
      assert(flowrowsigns[r] != 0);
      nnodes++;
      SCIP_CALL( addCommodityFlowRow(scip, flowrowsigns, flowrowcom, rowarcnodeid, capacityrowsigns, pluscom, minuscom,
                                     comcols, ncomcols, comcolssize, commerge,
                                     rows[r], flowrowsigns[r], 0, nnodes-1) );
      newrow = rows[r];
      
      /* big extraction loop */
      foundfit = TRUE;
      while( foundfit )
      {
         SCIP_COL** newrowcols;
         SCIP_Real* newrowvals;
         int newrowlen;
         int newrowscale;
         int nplussigns;
         int nminussigns;
         int nunknownplussigns;
         int nunknownminussigns;
         int i;

         assert(newrow != NULL);

         SCIP_CALL( printCommodities(scip, ncommodities, commerge, flowrowsigns, flowrowcom, colarcid, rowarcnodeid, pluscom, minuscom) );

         r = SCIProwGetLPPos(newrow);
         assert(r >= 0);
         newrowcols = SCIProwGetCols(newrow);
         newrowvals = SCIProwGetVals(newrow);
         newrowlen = SCIProwGetNNonz(newrow);

         /* 4. For all newly introduced arc flow variables of commodity k due to the new row:
          *    (a) find a capacity constraint which is compatible to all commodities:
          *         (i)   must have entry for arc flow variable
          *         (ii)  must not have entry for other arc flow variables in the commodity k
          *         (iii) must have at most (usually exactly) one entry for arc flow variables of each commodity
          *         (iv)  each entry that belongs to an arc flow variable of a commodity must not have the opposite
          *               sign in this commodity compared to commodity k
          *        If no capacity constraint has been found, goto 5.
          *    (b) For the capacity constraint, identify corresponding columns in other commodities:
          *         (i) For all potential arc flow variables in the capacity constraint that do not yet belong to a
          *             commodity:
          *              (1) For all flow conservation rows in which this column appears:
          *                  - Check if it fits to the new node flow conservation constraint in commodity k:
          *                     * columns already assigned to a commodity must have the same sign as in new node flow conservation constraint
          *                     * number of positive and negative coefficients of unassigned columns must match new node flow conservation constraint
          *                    Otherwise, ignore constraint.
          *                    Note that if the constraint has overlap with some commodity, it will automatically fit due to the design of the algorithm.
          *                  - If found flow conservation constraint has columns which are already assigned to some commodity,
          *                    add it to the same commodity and exit this loop.
          *                    Otherwise, remember row.
          *              (2) If no constraint was added to an existing commodity, generate new commodity and add the remembered row.
          */

         /* count number of plus and minus signs in new flow conservation constraint */
         if( (flowrowsigns[r] & LHSASSIGNED) != 0 )
            newrowscale = +1;
         else
         {
            assert((flowrowsigns[r] & RHSASSIGNED) != 0);
            newrowscale = -1;
         }
         nplussigns = 0;
         nminussigns = 0;
         nunknownplussigns = 0;
         nunknownminussigns = 0;
         for( i = 0; i < newrowlen; i++ )
         {
            int j;
            int arcid;

            j = SCIPcolGetLPPos(newrowcols[i]);
            if( j >= 0 )
               arcid = colarcid[j];
            else
               arcid = -1;
            assert(arcid < ncols);
            assert(arcid < narcs);

            if( newrowscale * newrowvals[i] > 0.0 )
            {
               nplussigns++;
               if( arcid >= 0 )
                  arcpattern[arcid] = +1;
               else
                  nunknownplussigns++;
            }
            else
            {
               nminussigns++;
               if( arcid >= 0 )
                  arcpattern[arcid] = -1;
               else
                  nunknownminussigns++;
            }
         }

         /* find capacity constraints for all new columns in the new row */
         for( i = 0; i < newrowlen; i++ )
         {
            c = SCIPcolGetLPPos(newrowcols[i]);
            if( c >= 0 && (pluscom[c] >= 0) + (minuscom[c] >= 0) == 1 ) /* exactly one of pluscom[c] and minuscom[c] is >= 0 */
            {
               SCIP_ROW** colrows;
               int collen;
               int j;

               SCIPdebugMessage("searching capacity row for column %d <%s>\n", c, SCIPvarGetName(SCIPcolGetVar(newrowcols[i])));

               /* find capacity constraint which is compatible to all commodities */
               colrows = SCIPcolGetRows(cols[c]);
               collen = SCIPcolGetNNonz(cols[c]);
               for( j = 0; j < collen; j++ )
               {
                  SCIP_ROW* row;
                  int capacityrowsign;
                  
                  row = colrows[j];
                  capacityrowsign = getCapacityRowFit(scip, ncommodities, capacityrowsigns, pluscom, minuscom, (pluscom[c] >= 0), row);
                  if( capacityrowsign != 0 )
                  {
                     /* update the arc pattern, because this column will be assigned to a new arc */
                     assert(arcpattern[narcs] == 0);
                     if( pluscom[c] >= 0 )
                     {
                        assert(nunknownplussigns >= 1);
                        arcpattern[narcs] = +1;
                        nunknownplussigns--;
                     }
                     else
                     {
                        assert(nunknownminussigns >= 1);
                        arcpattern[narcs] = -1;
                        nunknownminussigns--;
                     }

                     /* the capacity row defines a new arc */
                     narcs++;
                     assert(narcs <= ncols);

                     /* for the capacity constraint, identify corresponding columns in other commodities */
                     SCIP_CALL( addCapacityRow(scip, maxsigndistance, flowrowsigns, flowrowcom, colarcid, rowarcnodeid, capacityrowsigns, iscapacity,
                                               pluscom, minuscom, newrowscale,
                                               nplussigns, nminussigns, nunknownplussigns, nunknownminussigns, arcpattern,
                                               &comcols, &ncomcols, &comcolssize, &commerge, &ncommodities, &commoditiessize,
                                               c, row, capacityrowsign, nnodes-1, narcs-1) );
                     break;
                  }
               }
            }
         }

         /* reset arcpattern array */
         for( i = 0; i < newrowlen; i++ )
         {
            int j;
            int arcid;

            j = SCIPcolGetLPPos(newrowcols[i]);
            if( j < 0 )
               continue;

            arcid = colarcid[j];
            assert(arcid < ncols);
            assert(arcid < narcs);
            if( arcid >= 0 )
               arcpattern[arcid] = 0;
         }

         /* 6. For all commodities k: Find next node flow conservation constraint for commodity k.
          *    (a) For all columns already known in commodity k with only one entry:
          *        (i) For all flow conservation rows in which this column appears:
          *            - Check if it fits to this commodity. If yes, take it and goto 4.
          *    If nothing fits, exit.
          *
          * TODO: We do the extension of the commodities in a greedy fashion w.r.t. the order of the rows
          *       in the column vectors. It may be better to use some predefined priority order.
          */
         SCIPdebugMessage("searching next flow row\n");
         foundfit = FALSE;
         for( k = 0; k < ncommodities && !foundfit; k++ )
         {
            for( i = 0; i < ncomcols[k] && !foundfit; i++ )
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
                  for( j = 0; j < collen && !foundfit; j++ )
                  {
                     SCIP_ROW* row;
                     int flowrowsign;

                     row = colrows[j];

                     /* check if row fits into commodity */
                     flowrowsign = getFlowRowCommodityFit(flowrowsigns, flowrowcom, colarcid, iscapacity, pluscom, minuscom,
                                                          -1, -1, -1, -1, NULL, -1, +1, row, &k, NULL);
                     if( flowrowsign != 0 )
                     {
                        /* this flow conservation constraint defines a new node */
                        nnodes++;

                        /* row fits: insert into commodity */
                        SCIP_CALL( addCommodityFlowRow(scip, flowrowsigns, flowrowcom, rowarcnodeid, capacityrowsigns, pluscom, minuscom,
                                                       comcols, ncomcols, comcolssize, commerge, row, flowrowsign, k, nnodes-1) );
                        foundfit = TRUE;
                        newrow = row;
                     }
                  }
               }
            }
         }
      }

      flushCommodityMerges(nrows, commerge, flowrowcom);

#ifdef SCIP_DEBUG
      /* output commodities */
      SCIP_CALL( printCommodities(scip, ncommodities, commerge, flowrowsigns, flowrowcom, colarcid, rowarcnodeid, pluscom, minuscom) );

      /* output capacity constraints */
      printf("capacity constraints:\n");
      for( r = 0; r < nrows; r++ )
      {
         if( (capacityrowsigns[r] & (LHSASSIGNED | RHSASSIGNED)) != 0 )
            printf("%2d: <%s> [arc:%2d]\n", r, SCIProwGetName(rows[r]), rowarcnodeid[r]);
      }

      /* output unused constraints */
      printf("unused constraints:\n");
      for( r = 0; r < nrows; r++ )
      {
         if( (flowrowsigns[r] & (LHSASSIGNED | RHSASSIGNED)) == 0
             && (capacityrowsigns[r] & (LHSASSIGNED | RHSASSIGNED)) == 0 )
         {
            printf("%2d: <%s> [flowsign=%d capsign=%d]\n", r, SCIProwGetName(rows[r]), flowrowsigns[r], capacityrowsigns[r]);
            /*SCIP_CALL( SCIPprintRow(scip, rows[r], NULL) );*/
         }
      }
      printf("\n");

      /* output unused columns */
      printf("unused columns:\n");
      for( c = 0; c < ncols; c++ )
      {
         if( pluscom[c] == -1 && minuscom[c] == -1 )
         {
            SCIP_VAR* var = SCIPcolGetVar(cols[c]);
            printf("<%s> [%g,%g]\n", SCIPvarGetName(var), SCIPvarGetLbGlobal(var), SCIPvarGetUbGlobal(var));
         }
      }
#endif

      /* free memory */
      for( k = 0; k < ncommodities; k++ )
         SCIPfreeMemoryArray(scip, &comcols[k]);
      SCIPfreeMemoryArray(scip, &comcols);
      SCIPfreeMemoryArray(scip, &ncomcols);
      SCIPfreeMemoryArray(scip, &comcolssize);
      SCIPfreeMemoryArray(scip, &commerge);

      SCIPfreeBufferArray(scip, &arcpattern);
      SCIPfreeBufferArray(scip, &rowarcnodeid);
      SCIPfreeBufferArray(scip, &colarcid);
      SCIPfreeBufferArray(scip, &flowrowcom);
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
      SCIP_CALL( mcfnetworkExtract(scip, sepadata->maxsigndistance, &sepadata->mcfnetwork) );
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
   SCIP_CALL( SCIPaddIntParam(scip,
         "separating/mcf/maxsigndistance",
         "maximum Hamming distance of flow conservation constraint sign patterns of the same node",
         &sepadata->maxsigndistance, TRUE, DEFAULT_MAXSIGNDISTANCE, 0, INT_MAX, NULL, NULL) );

   return SCIP_OKAY;
}
