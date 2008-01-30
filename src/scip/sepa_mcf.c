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
#pragma ident "@(#) $Id: sepa_mcf.c,v 1.12 2008/01/30 17:17:05 bzfpfend Exp $"

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

/** internal MCF extraction data to pass to subroutines */
struct mcfdata
{
   unsigned char*        flowrowsigns;       /**< potential or actual sides of rows to be used as flow conservation constraint */
   SCIP_Real*            flowrowscalars;     /**< scalar of rows to transform into +/-1 coefficients */
   SCIP_Real*            flowrowscores;      /**< score value indicating how sure we are that this is indeed a flow conservation constraint */
   unsigned char*        capacityrowsigns;   /**< potential or actual sides of rows to be used as capacity constraint */
   SCIP_Real*            capacityrowscores;  /**< score value indicating how sure we are that this is indeed a capacity constraint */
   int*                  flowcands;          /**< list of row indices that are candidates for flow conservation constraints */
   int                   nflowcands;         /**< number of elements in flow candidate list */
   int*                  capacitycands;      /**< list of row indices that are candidates for capacity constraints */
   int                   ncapacitycands;     /**< number of elements in capacity candidate list */
   SCIP_Bool*            plusflow;           /**< is column c member of a flow row with coefficient +1? */
   SCIP_Bool*            minusflow;          /**< is column c member of a flow row with coefficient -1? */
   int                   ncommodities;       /**< number of commodities */
   int*                  commoditysigns;     /**< +1: regular, -1: all arcs have opposite direction; 0: undecided */
   int                   commoditysignssize; /**< size of commoditysigns array */
   int*                  colcommodity;       /**< commodity number of each column, or -1 */
   int*                  rowcommodity;       /**< commodity number of each row, or -1 */
   int*                  colarcid;           /**< arc id of each flow column, or -1 */
   int*                  rowarcid;           /**< arc id of each capacity row, or -1 */
   int*                  rownodeid;          /**< node id of each flow conservation row, or -1 */
   int*                  newcols;            /**< columns of current commodity that have to be inspected for incident flow conservation rows */
   int                   nnewcols;           /**< number of newcols */
   int                   narcs;              /**< number of arcs in the extracted graph */
   int                   nnodes;             /**< number of nodes in the extracted graph */
   SCIP_ROW**            capacityrows;       /**< capacity row for each arc */
   int                   capacityrowssize;   /**< size of array */
   SCIP_Bool*            colisincident;      /**< temporary memory for column collection */
};
typedef struct mcfdata MCFDATA;              /**< internal MCF extraction data to pass to subroutines */




/*
 * Local methods
 */

#define LHSPOSSIBLE  1
#define RHSPOSSIBLE  2
#define LHSASSIGNED  4
#define RHSASSIGNED  8
#define DISCARDED   16

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

#ifdef SCIP_DEBUG
/** displays commodities and its members */
static
void printCommodities(
   SCIP*                 scip,               /**< SCIP data structure */ 
   MCFDATA*              mcfdata             /**< internal MCF extraction data to pass to subroutines */
   )
{
   unsigned char* capacityrowsigns = mcfdata->capacityrowsigns;
   int            ncommodities     = mcfdata->ncommodities;
   int*           commoditysigns   = mcfdata->commoditysigns;
   int*           colcommodity     = mcfdata->colcommodity;
   int*           rowcommodity     = mcfdata->rowcommodity;
   int*           colarcid         = mcfdata->colarcid;
   int*           rownodeid        = mcfdata->rownodeid;
   SCIP_ROW**     capacityrows     = mcfdata->capacityrows;

   SCIP_COL** cols;
   SCIP_ROW** rows;
   int ncols;
   int nrows;
   int k;
   int c;
   int r;
   int a;

   cols = SCIPgetLPCols(scip);
   ncols = SCIPgetNLPCols(scip);
   rows = SCIPgetLPRows(scip);
   nrows = SCIPgetNLPRows(scip);

   for( k = 0; k < ncommodities; k++ )
   {
      printf("commodity %d (sign: %+d):\n", k, commoditysigns[k]);

      for( c = 0; c < ncols; c++ )
      {
         if( colcommodity[c] == k )
            printf(" col <%s>: arc %d\n", SCIPvarGetName(SCIPcolGetVar(cols[c])), colarcid != NULL ? colarcid[c] : -1);
      }
      for( r = 0; r < nrows; r++ )
      {
         if( rowcommodity[r] == k )
            printf(" row <%s>: node %d\n", SCIProwGetName(rows[r]), rownodeid != NULL ? rownodeid[r] : -1);
      }
      printf("\n");
   }

   printf("capacities:\n");
   for( a = 0; a < mcfdata->narcs; a++ )
   {
      printf("  arc %d: ", a);
      if( capacityrows[a] != NULL )
      {
         r = SCIProwGetLPPos(capacityrows[a]);
         assert(0 <= r && r < nrows);
         if( (capacityrowsigns[r] & LHSASSIGNED) != 0 )
            printf(" row <%s> [sign:-1]\n", SCIProwGetName(rows[r]));
         else if( (capacityrowsigns[r] & RHSASSIGNED) != 0 )
            printf(" row <%s> [sign:+1]\n", SCIProwGetName(rows[r]));
      }
      else
         printf(" -\n");
   }
   printf("\n");

   printf("unused columns:\n");
   for( c = 0; c < ncols; c++ )
   {
      if( colcommodity[c] == -1 )
      {
         SCIP_VAR* var = SCIPcolGetVar(cols[c]);
         printf(" col <%s> [%g,%g]\n", SCIPvarGetName(var), SCIPvarGetLbGlobal(var), SCIPvarGetUbGlobal(var));
      }
   }
   printf("\n");

   printf("unused rows:\n");
   for( r = 0; r < nrows; r++ )
   {
      if( rowcommodity[r] == -1 && (capacityrowsigns[r] & (LHSASSIGNED | RHSASSIGNED)) == 0 )
         printf(" row <%s>\n", SCIProwGetName(rows[r]));
   }
   printf("\n");
}
#endif

/** comparator method for flow and capacity row candidates */
static
SCIP_DECL_SORTINDCOMP(compCands)
{
   SCIP_Real* rowscores = (SCIP_Real*)dataptr;

   return rowscores[ind2] - rowscores[ind1];
}

/** extracts flow conservation and capacity rows from the LP */
static
SCIP_RETCODE extractRows(
   SCIP*                 scip,               /**< SCIP data structure */ 
   MCFDATA*              mcfdata             /**< internal MCF extraction data to pass to subroutines */
   )
{
   unsigned char* flowrowsigns;
   SCIP_Real*     flowrowscalars;
   SCIP_Real*     flowrowscores;
   unsigned char* capacityrowsigns;
   SCIP_Real*     capacityrowscores;
   int*           flowcands;
   int*           capacitycands;

   SCIP_ROW** rows;
   int nrows;
   int r;

   SCIP_CALL( SCIPgetLPRowsData(scip, &rows, &nrows) );

   /* allocate temporary memory for extraction data */
   SCIP_CALL( SCIPallocMemoryArray(scip, &mcfdata->flowrowsigns, nrows) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &mcfdata->flowrowscalars, nrows) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &mcfdata->flowrowscores, nrows) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &mcfdata->capacityrowsigns, nrows) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &mcfdata->capacityrowscores, nrows) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &mcfdata->flowcands, nrows) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &mcfdata->capacitycands, nrows) );
   flowrowsigns      = mcfdata->flowrowsigns;
   flowrowscalars    = mcfdata->flowrowscalars;
   flowrowscores     = mcfdata->flowrowscores;
   capacityrowsigns  = mcfdata->capacityrowsigns;
   capacityrowscores = mcfdata->capacityrowscores;
   flowcands         = mcfdata->flowcands;
   capacitycands     = mcfdata->capacitycands;

   for( r = 0; r < nrows; r++ )
   {
      SCIP_ROW* row;
      SCIP_COL** rowcols;
      SCIP_Real* rowvals;
      SCIP_Real rowlhs;
      SCIP_Real rowrhs;
      int rowlen;
      int nbinvars;
      int nintvars;
      int nimplintvars;
      int ncontvars;
      SCIP_Real coef;
      SCIP_Bool hasposcoef;
      SCIP_Bool hasnegcoef;
      SCIP_Bool hasposcontcoef;
      SCIP_Bool hasnegcontcoef;
      SCIP_Bool hasposintcoef;
      SCIP_Bool hasnegintcoef;
      SCIP_Bool hasinteger;
      SCIP_Real samecontcoef;
      unsigned int rowsign;
      int i;

      row = rows[r];
      assert(SCIProwGetLPPos(row) == r);

      flowrowsigns[r] = 0;
      flowrowscalars[r] = 0.0;
      flowrowscores[r] = 0.0;
      capacityrowsigns[r] = 0;
      capacityrowscores[r] = 0.0;

      rowlen = SCIProwGetNNonz(row);
      if( rowlen == 0 )
         continue;
      rowcols = SCIProwGetCols(row);
      rowvals = SCIProwGetVals(row);
      rowlhs = SCIProwGetLhs(row);
      rowrhs = SCIProwGetRhs(row);

      /* identify flow conservation constraints */
      coef = ABS(rowvals[0]);
      hasposcoef = FALSE;
      hasnegcoef = FALSE;
      nbinvars = 0;
      nintvars = 0;
      nimplintvars = 0;
      ncontvars = 0;
      for( i = 0; i < rowlen; i++ )
      {
         SCIP_Real absval = ABS(rowvals[i]);
         if( !SCIPisEQ(scip, absval, coef) )
            break;

#if 0
         if( SCIPvarGetType(SCIPcolGetVar(rowcols[i])) != SCIP_VARTYPE_CONTINUOUS )
            break; /*????????????????????*/
#endif
         hasposcoef = hasposcoef || (rowvals[i] > 0.0);
         hasnegcoef = hasnegcoef || (rowvals[i] < 0.0);
         switch( SCIPvarGetType(SCIPcolGetVar(rowcols[i])) )
         {
         case SCIP_VARTYPE_BINARY:
            nbinvars++;
            break;
         case SCIP_VARTYPE_INTEGER:
            nintvars++;
            break;
         case SCIP_VARTYPE_IMPLINT:
            nimplintvars++;
            break;
         case SCIP_VARTYPE_CONTINUOUS:
            ncontvars++;
            break;
         default:
            SCIPerrorMessage("unknown variable type\n");
            SCIPABORT();
         }
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
         flowcands[mcfdata->nflowcands] = r;
         mcfdata->nflowcands++;
      }

      /* identify capacity constraints */
      hasposcontcoef = FALSE;
      hasnegcontcoef = FALSE;
      hasposintcoef = FALSE;
      hasnegintcoef = FALSE;
      hasinteger = FALSE;
      samecontcoef = 0.0;
      rowsign = 0;
      if( !SCIPisInfinity(scip, -rowlhs) )
         rowsign |= LHSPOSSIBLE;
      if( !SCIPisInfinity(scip, rowrhs) )
         rowsign |= RHSPOSSIBLE;
      for( i = 0; i < rowlen && rowsign != 0; i++ )
      {
         if( SCIPvarGetType(SCIPcolGetVar(rowcols[i])) == SCIP_VARTYPE_CONTINUOUS )
         {
            if( samecontcoef == 0.0 )
               samecontcoef = rowvals[i];
            else if( !SCIPisEQ(scip, samecontcoef, rowvals[i]) )
               samecontcoef = SCIP_REAL_MAX;
            
            if( rowvals[i] > 0.0 )
            {
               rowsign &= ~LHSPOSSIBLE;
               hasposcontcoef = TRUE;
            }
            else
            {
               rowsign &= ~RHSPOSSIBLE;
               hasnegcontcoef = TRUE;
            }
         }
         else
         {
            hasinteger = TRUE;
            if( rowvals[i] > 0.0 )
               hasposintcoef = TRUE;
            else
               hasnegintcoef = TRUE;
         }
      }
      if( i == rowlen && hasinteger && rowsign != 0 )
      {
         capacityrowsigns[r] = rowsign;
         capacitycands[mcfdata->ncapacitycands] = r;
         mcfdata->ncapacitycands++;
      }

      /* calculate flow row score */
      if( (flowrowsigns[r] & (LHSPOSSIBLE | RHSPOSSIBLE)) != 0 )
      {
         /* row is an equation: score +10 */
         if( (flowrowsigns[r] & (LHSPOSSIBLE | RHSPOSSIBLE)) == (LHSPOSSIBLE | RHSPOSSIBLE) )
            flowrowscores[r] += 10.0;

         /* row is not a capacity constraint: score +10 */
         if( (capacityrowsigns[r] & (LHSPOSSIBLE | RHSPOSSIBLE)) == 0 )
            flowrowscores[r] += 10.0;

         /* row does not need to be scaled: score +10 */
         if( SCIPisEQ(scip, flowrowscalars[r], 1.0) )
            flowrowscores[r] += 10.0;

         /* row has positive and negative coefficients: score +5 */
         if( hasposcoef && hasnegcoef )
            flowrowscores[r] += 5.0;

         /* all variables are of the same type:
          *    continuous: score +10
          *    integer:    score  +5
          *    binary:     score  +1
          */
         if( ncontvars == rowlen )
            flowrowscores[r] += 10.0;
         else if( nintvars + nimplintvars == rowlen )
            flowrowscores[r] += 5.0;
         else if( nbinvars == rowlen )
            flowrowscores[r] += 1.0;

         /**@todo go through list of several model types, depending on the current model type throw away invalid constraints
          *       instead of assigning a low score
          */
      }

      /* calculate capacity row score */
      if( (capacityrowsigns[r] & (LHSPOSSIBLE | RHSPOSSIBLE)) != 0 )
      {
         /* row is of type f - c*x <= b: score +10 */
         if( (capacityrowsigns[r] & RHSPOSSIBLE) != 0 && hasposcontcoef && !hasnegcontcoef && !hasposintcoef && hasnegintcoef )
            capacityrowscores[r] += 10.0;
         if( (capacityrowsigns[r] & LHSPOSSIBLE) != 0 && !hasposcontcoef && hasnegcontcoef && hasposintcoef && !hasnegintcoef )
            capacityrowscores[r] += 10.0;

         /* all coefficients of continuous variables are +1 or all are -1: score +5 */
         if( SCIPisEQ(scip, ABS(samecontcoef), 1.0) )
            capacityrowscores[r] += 5.0;

         /* all coefficients of continuous variables are equal: score +2 */
         if( samecontcoef != 0.0 && samecontcoef != SCIP_REAL_MAX )
            capacityrowscores[r] += 2.0;

         /* row is a <= row with non-negative right hand side: score +1 */
         if( (capacityrowsigns[r] & RHSPOSSIBLE) != 0 && !SCIPisNegative(scip, rowrhs)  )
            capacityrowscores[r] += 1.0;

         /* row is an inequality: score +1 */
         if( SCIPisInfinity(scip, -rowlhs) != SCIPisInfinity(scip, rowrhs) )
            capacityrowscores[r] += 1.0;
      }
   }

   /* sort candidates by score */
   SCIPbsortInd((void*)flowrowscores, mcfdata->flowcands, mcfdata->nflowcands, compCands);
   SCIPbsortInd((void*)capacityrowscores, mcfdata->capacitycands, mcfdata->ncapacitycands, compCands);

   SCIPdebugMessage("flow conservation candidates:\n");
   for( r = 0; r < mcfdata->nflowcands; r++ )
   {
      //SCIPdebug(SCIPprintRow(scip, rows[mcfdata->flowcands[r]], NULL));
      SCIPdebugMessage("%4d [score: %2g]: %s\n", r, flowrowscores[mcfdata->flowcands[r]], SCIProwGetName(rows[mcfdata->flowcands[r]]));
   }
   SCIPdebugMessage("capacity candidates:\n");
   for( r = 0; r < mcfdata->ncapacitycands; r++ )
   {
      //SCIPdebug(SCIPprintRow(scip, rows[mcfdata->capacitycands[r]], NULL));
      SCIPdebugMessage("%4d [score: %2g]: %s\n", r, capacityrowscores[mcfdata->capacitycands[r]], SCIProwGetName(rows[mcfdata->capacitycands[r]]));
   }

   return SCIP_OKAY;
}

/** creates a new commodity */
static
SCIP_RETCODE createNewCommodity(
   SCIP*                 scip,               /**< SCIP data structure */ 
   MCFDATA*              mcfdata             /**< internal MCF extraction data to pass to subroutines */
   )
{
   /* get memory for commoditysigns array */
   assert(mcfdata->ncommodities <= mcfdata->commoditysignssize);
   if( mcfdata->ncommodities == mcfdata->commoditysignssize )
   {
      mcfdata->commoditysignssize = MAX(2*mcfdata->commoditysignssize, mcfdata->ncommodities+1);
      SCIP_CALL( SCIPreallocMemoryArray(scip, &mcfdata->commoditysigns, mcfdata->commoditysignssize) );
   }
   assert(mcfdata->ncommodities < mcfdata->commoditysignssize);

   /* create commodity */
   SCIPdebugMessage("**** creating new commodity %d ****\n", mcfdata->ncommodities);
   mcfdata->commoditysigns[mcfdata->ncommodities] = 0;
   mcfdata->ncommodities++;

   return SCIP_OKAY;
}

/** adds the given flow row and all involved columns to the current commodity */
static
void addFlowrowToCommodity(
   SCIP*                 scip,               /**< SCIP data structure */ 
   MCFDATA*              mcfdata,            /**< internal MCF extraction data to pass to subroutines */
   SCIP_ROW*             row,                /**< flow row to add to current commodity */
   unsigned char         rowsign             /**< possible flow row signs to use */
   )
{
   unsigned char* flowrowsigns   = mcfdata->flowrowsigns;
   SCIP_Bool*     plusflow       = mcfdata->plusflow;
   SCIP_Bool*     minusflow      = mcfdata->minusflow;
   int            ncommodities   = mcfdata->ncommodities;
   int*           commoditysigns = mcfdata->commoditysigns;
   int*           colcommodity   = mcfdata->colcommodity;
   int*           rowcommodity   = mcfdata->rowcommodity;
   int*           newcols        = mcfdata->newcols;

   SCIP_COL** rowcols;
   SCIP_Real* rowvals;
   int rowlen;
   int rowscale;
   int r;
   int k;
   int i;

   k = ncommodities-1;
   assert(k >= 0);

   r = SCIProwGetLPPos(row);
   assert(r >= 0);
   assert(rowcommodity[r] == -1);
   assert((flowrowsigns[r] | rowsign) == flowrowsigns[r]);
   assert((rowsign & (LHSPOSSIBLE | RHSPOSSIBLE)) == rowsign);
   assert(rowsign != 0);

   /* decide the sign (direction) of the row */
   if( rowsign == LHSPOSSIBLE )
   {
      rowsign = LHSASSIGNED;
      commoditysigns[k] = +1; /* we cannot switch directions */
   }
   else if( rowsign == RHSPOSSIBLE )
   {
      rowsign = RHSASSIGNED;
      commoditysigns[k] = +1; /* we cannot switch directions */
   }
   else
   {
      assert(rowsign == (LHSPOSSIBLE | RHSPOSSIBLE));
      rowsign = LHSASSIGNED; /* for now, choose lhs */
   }
   if( rowsign == LHSASSIGNED )
      rowscale = +1;
   else
      rowscale = -1;
   flowrowsigns[r] |= rowsign;

   SCIPdebugMessage("adding flow row %d <%s> with sign %d to commodity %d\n", r, SCIProwGetName(row), rowscale, k);
   SCIPdebug( SCIPprintRow(scip, row, NULL); );

   /* add row to commodity */
   rowcommodity[r] = k;
   rowcols = SCIProwGetCols(row);
   rowvals = SCIProwGetVals(row);
   rowlen = SCIProwGetNLPNonz(row);
   for( i = 0; i < rowlen; i++ )
   {
      SCIP_Real val;
      int c;

      c = SCIPcolGetLPPos(rowcols[i]);
      assert(0 <= c && c < SCIPgetNLPCols(scip));

      /* assign column to commodity */
      if( colcommodity[c] == -1 )
      {
         assert(!plusflow[c]);
         assert(!minusflow[c]);
         assert(mcfdata->nnewcols < SCIPgetNLPCols(scip));
         colcommodity[c] = k;
         SCIPdebugMessage(" -> push x%d to stack at pos %d\n", c, mcfdata->nnewcols); /*??????????????????*/
         newcols[mcfdata->nnewcols] = c;
         mcfdata->nnewcols++;
      }
      assert(colcommodity[c] == k);

      /* update plusflow/minusflow */
      val = rowscale * rowvals[i];
      if( val > 0.0 )
      {
         assert(!plusflow[c]);
         plusflow[c] = TRUE;
      }
      else
      {
         assert(!minusflow[c]);
         minusflow[c] = TRUE;
      }
   }
}

/** checks whether the given row fits into the given commodity and returns the possible flow row signs */
static
void getFlowrowFit(
   SCIP*                 scip,               /**< SCIP data structure */ 
   MCFDATA*              mcfdata,            /**< internal MCF extraction data to pass to subroutines */
   SCIP_ROW*             row,                /**< flow row to check */
   int                   k,                  /**< commodity that the flow row should enter */
   unsigned char*        rowsign             /**< pointer to store the possible flow row signs */
   )
{
   unsigned char* flowrowsigns = mcfdata->flowrowsigns;
   SCIP_Bool*     plusflow     = mcfdata->plusflow;
   SCIP_Bool*     minusflow    = mcfdata->minusflow;
   int*           colcommodity = mcfdata->colcommodity;
   int*           rowcommodity = mcfdata->rowcommodity;

   SCIP_COL** rowcols;
   SCIP_Real* rowvals;
   int rowlen;
   int flowrowsign;
   int r;
   int j;

   *rowsign = 0;

   r = SCIProwGetLPPos(row);
   assert(0 <= r && r < SCIPgetNLPRows(scip));

   /* ignore rows that are already used */
   if( rowcommodity[r] != -1 )
      return;

   /* check if row is an available flow row */
   flowrowsign = flowrowsigns[r];
   assert((flowrowsign & (LHSPOSSIBLE | RHSPOSSIBLE | DISCARDED)) == flowrowsign);
   if( (flowrowsign & DISCARDED) != 0 )
      return;
   if( (flowrowsign & (LHSPOSSIBLE | RHSPOSSIBLE)) == 0 )
      return;

   /* check whether the row fits w.r.t. the signs of the coefficients */
   rowcols = SCIProwGetCols(row);
   rowvals = SCIProwGetVals(row);
   rowlen = SCIProwGetNLPNonz(row);
   for( j = 0; j < rowlen && flowrowsign != 0; j++ )
   {
      int rowc;

      rowc = SCIPcolGetLPPos(rowcols[j]);
      assert(0 <= rowc && rowc < SCIPgetNLPCols(scip));

      /* check if column already belongs to the same commodity */
      if( colcommodity[rowc] == k )
      {
         /* column only fits if it is not yet present with the same sign */
         if( plusflow[rowc] )
         {
            /* column must not be included with positive sign */
            if( rowvals[j] > 0.0 )
               flowrowsign &= ~LHSPOSSIBLE;
            else
               flowrowsign &= ~RHSPOSSIBLE;
         }
         if( minusflow[rowc] )
         {
            /* column must not be included with negative sign */
            if( rowvals[j] > 0.0 )
               flowrowsign &= ~RHSPOSSIBLE;
            else
               flowrowsign &= ~LHSPOSSIBLE;
         }
      }
      else if( colcommodity[rowc] != -1 )
      {
         /* column does not fit if it already belongs to a different commodity */
         flowrowsign = 0;
      }
   }

   if( flowrowsign == 0 )
   {
      /* we can discard the row, since it can also not be member of a different commodity */
      SCIPdebugMessage(" -> discard flow row %d <%s>\n", r, SCIProwGetName(row));
      flowrowsigns[r] |= DISCARDED;
   }

   *rowsign = flowrowsign;
}

/** returns a flow conservation row that fits into the current commodity, or NULL */
static
void getNextFlowrow(
   SCIP*                 scip,               /**< SCIP data structure */ 
   MCFDATA*              mcfdata,            /**< internal MCF extraction data to pass to subroutines */
   SCIP_ROW**            nextrow,            /**< pointer to store next row */
   unsigned char*        nextrowsign         /**< pointer to store possible signs of next row */
   )
{
   SCIP_Bool* plusflow     = mcfdata->plusflow;
   SCIP_Bool* minusflow    = mcfdata->minusflow;
   int*       newcols      = mcfdata->newcols;
   int        ncommodities = mcfdata->ncommodities;

   SCIP_COL** cols;
   int k;

   assert(nextrow != NULL);
   assert(nextrowsign != NULL);

   *nextrow = NULL;
   *nextrowsign = 0;

   k = ncommodities-1;

   cols = SCIPgetLPCols(scip);
   assert(cols != NULL);

   /* check if there are any columns left in the commodity that have not yet been inspected for incident flow rows */
   while( mcfdata->nnewcols > 0 && *nextrow == NULL )
   {
      SCIP_COL* col;
      SCIP_ROW** colrows;
      int collen;
      int c;
      int i;

      /* pop next new column from stack */
      c = newcols[mcfdata->nnewcols-1];
      mcfdata->nnewcols--;
      SCIPdebugMessage(" -> pop x%d from stack at pos %d\n", c, mcfdata->nnewcols); /*??????????????????*/
      assert(0 <= c && c < SCIPgetNLPCols(scip));

      /* check if this columns already as both signs */
      assert(plusflow[c] || minusflow[c]);
      if( plusflow[c] && minusflow[c] )
         continue;

      /* check whether column is incident to a valid flow row that fits into the current commodity */
      col = cols[c];
      colrows = SCIPcolGetRows(col);
      collen = SCIPcolGetNLPNonz(col);
      for( i = 0; i < collen; i++ )
      {
         SCIP_ROW* row;
         unsigned char flowrowsign;

         row = colrows[i];

         /* check if row fits into the current commodity */
         getFlowrowFit(scip, mcfdata, row, k, &flowrowsign);

         /* do we have a winner? */
         if( flowrowsign != 0 )
         {
            *nextrow = row;
            *nextrowsign = flowrowsign;
            break;
         }
      }
   }
}

/** extracts flow conservation rows and puts them into commodities */
static
SCIP_RETCODE extractFlow(
   SCIP*                 scip,               /**< SCIP data structure */ 
   MCFDATA*              mcfdata             /**< internal MCF extraction data to pass to subroutines */
   )
{
   int* flowcands = mcfdata->flowcands;

   SCIP_Bool* plusflow;
   SCIP_Bool* minusflow;
   int* colcommodity;
   int* rowcommodity;

   SCIP_ROW** rows;
   int nrows;
   int ncols;
   int i;
   int c;
   int r;

   /* get LP data */
   rows = SCIPgetLPRows(scip);
   nrows = SCIPgetNLPRows(scip);
   ncols = SCIPgetNLPCols(scip);

   /* allocate temporary memory */
   SCIP_CALL( SCIPallocMemoryArray(scip, &mcfdata->plusflow, ncols) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &mcfdata->minusflow, ncols) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &mcfdata->colcommodity, ncols) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &mcfdata->rowcommodity, nrows) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &mcfdata->newcols, ncols) );
   plusflow = mcfdata->plusflow;
   minusflow = mcfdata->minusflow;
   colcommodity = mcfdata->colcommodity;
   rowcommodity = mcfdata->rowcommodity;

   /* 3. Extract network structure of flow conservation constraints:
    *    (a) Initialize plusflow[c] = minusflow[c] = FALSE for all columns c and other local data.
    */
   BMSclearMemoryArray(plusflow, ncols);
   BMSclearMemoryArray(minusflow, ncols);
   for( c = 0; c < ncols; c++ )
      colcommodity[c] = -1;
   for( r = 0; r < nrows; r++ )
      rowcommodity[r] = -1;

   /*    (b) As long as there are flow conservation candidates left:
    *        (i) Create new commodity and use first flow conservation constraint as newrow.
    *       (ii) Add newrow to commodity, update pluscom/minuscom accordingly.
    *      (iii) For the newly added columns search for an incident flow conservation constraint. Pick the one of highest ranking.
    *       (iv) If found, set newrow to this row and goto (ii).
    */
   for( i = 0; i < mcfdata->nflowcands; i++ )
   {
      SCIP_ROW* newrow;
      unsigned char newrowsign;

      r = flowcands[i];
      assert(0 <= r && r < nrows);
      newrow = rows[r];

      /* check if row fits into a new commodity */
      getFlowrowFit(scip, mcfdata, newrow, mcfdata->ncommodities, &newrowsign);
      if( newrowsign == 0 )
         continue;

      /* start new commodity */
      SCIP_CALL( createNewCommodity(scip, mcfdata) );

      /* fill commodity with flow conservation constraints */
      do
      {
         /* add new row to commodity */
         addFlowrowToCommodity(scip, mcfdata, newrow, newrowsign);

         /* get next row to add */
         getNextFlowrow(scip, mcfdata, &newrow, &newrowsign);
      }
      while( newrow != NULL );
   }

   return SCIP_OKAY;
}

/** identifies capacity constraints for the arcs and assigns arc ids to columns and capacity constraints */
static
SCIP_RETCODE extractCapacities(
   SCIP*                 scip,               /**< SCIP data structure */ 
   MCFDATA*              mcfdata             /**< internal MCF extraction data to pass to subroutines */
   )
{
   unsigned char* capacityrowsigns  = mcfdata->capacityrowsigns;
   SCIP_Real*     capacityrowscores = mcfdata->capacityrowscores;
   int*           colcommodity      = mcfdata->colcommodity;
   int*           rowcommodity      = mcfdata->rowcommodity;

   int* colarcid;
   int* rowarcid;

   SCIP_ROW** rows;
   SCIP_COL** cols;
   int nrows;
   int ncols;
   int r;
   int c;

   assert(mcfdata->narcs == 0);

   /* get LP data */
   SCIP_CALL( SCIPgetLPRowsData(scip, &rows, &nrows) );
   SCIP_CALL( SCIPgetLPColsData(scip, &cols, &ncols) );

   /* allocate temporary memory for extraction data */
   SCIP_CALL( SCIPallocMemoryArray(scip, &mcfdata->colarcid, ncols) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &mcfdata->rowarcid, nrows) );
   colarcid = mcfdata->colarcid;
   rowarcid = mcfdata->rowarcid;

   /* initialize arcid arrays */
   for( c = 0; c < ncols; c++ )
      colarcid[c] = -1;
   for( r = 0; r < nrows; r++ )
      rowarcid[r] = -1;

   /* for each column, search for a capacity constraint */
   for( c = 0; c < ncols; c++ )
   {
      SCIP_ROW* bestcapacityrow;
      SCIP_Real bestscore;
      SCIP_ROW** colrows;
      int collen;
      int i;

      /* ignore columns that are not flow variables */
      if( colcommodity[c] == -1 )
         continue;

      /* ignore columns that are already assigned to an arc */
      if( colarcid[c] >= 0 )
         continue;

      /* scan the column to search for valid capacity constraints */
      bestcapacityrow = NULL;
      bestscore = 0.0;
      colrows = SCIPcolGetRows(cols[c]);
      collen = SCIPcolGetNLPNonz(cols[c]);
      for( i = 0; i < collen; i++ )
      {
         r = SCIProwGetLPPos(colrows[i]);
         assert(0 <= r && r < nrows);

         /* row must not be already assigned */
         assert((capacityrowsigns[r] & (LHSPOSSIBLE | RHSPOSSIBLE)) == capacityrowsigns[r]);

         /* ignore rows that are not capacity candidates */
         if( (capacityrowsigns[r] & (LHSPOSSIBLE | RHSPOSSIBLE)) == 0 )
            continue;

         /* ignore rows that are already used as flow conservation constraints */
         if( rowcommodity[r] != -1 )
            continue;

         /* check if this capacity candidate has better score */
         if( capacityrowscores[r] > bestscore )
         {
            bestcapacityrow = colrows[i];
            bestscore = capacityrowscores[r];
         }
      }

      /* if no capacity row has been found, leave the column unassigned */
      if( bestcapacityrow != NULL )
      {
         SCIP_COL** rowcols;
         int rowlen;

         /* store the row */
         assert(mcfdata->narcs <= mcfdata->capacityrowssize);
         if( mcfdata->narcs == mcfdata->capacityrowssize )
         {
            mcfdata->capacityrowssize = MAX(2*mcfdata->capacityrowssize, mcfdata->narcs+1);
            SCIP_CALL( SCIPreallocMemoryArray(scip, &mcfdata->capacityrows, mcfdata->capacityrowssize) );
         }
         assert(mcfdata->narcs < mcfdata->capacityrowssize);
         mcfdata->capacityrows[mcfdata->narcs] = bestcapacityrow;

         /* assign the capacity row to a new arc id */
         r = SCIProwGetLPPos(bestcapacityrow);
         assert(0 <= r && r < nrows);
         rowarcid[r] = mcfdata->narcs;

         /* decide which sign to use */
         if( (capacityrowsigns[r] & RHSPOSSIBLE) != 0 )
            capacityrowsigns[r] |= RHSASSIGNED;
         else
         {
            assert((capacityrowsigns[r] & LHSPOSSIBLE) != 0);
            capacityrowsigns[r] |= LHSASSIGNED;
         }

         SCIPdebugMessage("assigning capacity row %d <%s> with sign %+d to arc %d\n",
                          r, SCIProwGetName(bestcapacityrow), (capacityrowsigns[r] & RHSASSIGNED) != 0 ? +1 : -1, mcfdata->narcs);

         /* assign all involved flow variables to the new arc id */
         SCIPdebugMessage(" -> flow:");
         rowcols = SCIProwGetCols(bestcapacityrow);
         rowlen = SCIProwGetNLPNonz(bestcapacityrow);
         for( i = 0; i < rowlen; i++ )
         {
            int rowc;

            rowc = SCIPcolGetLPPos(rowcols[i]);
            assert(0 <= rowc && rowc < ncols);

            if( colcommodity[rowc] >= 0 )
            {
               SCIPdebug( printf(" x%d<%s>[%d]", rowc, SCIPvarGetName(SCIPcolGetVar(rowcols[i])), colcommodity[rowc]) );
               colarcid[rowc] = mcfdata->narcs;
            }
         }
         SCIPdebug( printf("\n") );

         /* increase number of arcs */
         mcfdata->narcs++;
      }
      else
      {
         SCIPdebugMessage("no capacity row found for column x%d <%s> in commodity %d\n", c, SCIPvarGetName(SCIPcolGetVar(cols[c])), colcommodity[c]);
      }
   }

   return SCIP_OKAY;
}

/** collects all flow columns of all commodities (except the one of the base row) that are incident to the node described by the given flow row */
static
void collectIncidentFlowCols(
   SCIP*                 scip,               /**< SCIP data structure */ 
   MCFDATA*              mcfdata,            /**< internal MCF extraction data to pass to subroutines */
   SCIP_ROW*             flowrow,            /**< flow conservation constraint that defines the node */
   int                   basecommodity       /**< commodity of the base row */
   )
{
   int*           colcommodity  = mcfdata->colcommodity;
   int*           colarcid      = mcfdata->colarcid;
   int*           newcols       = mcfdata->newcols;
   SCIP_ROW**     capacityrows  = mcfdata->capacityrows;
   SCIP_Bool*     colisincident = mcfdata->colisincident;

   SCIP_COL** rowcols;
   int rowlen;
   int i;

#ifndef NDEBUG
   /* check that the marker array is correctly initialized */
   for( i = 0; i < SCIPgetNLPCols(scip); i++ )
      assert(!colisincident[i]);
#endif

   /* loop through all flow columns in the flow conservation constraint */
   rowcols = SCIProwGetCols(flowrow);
   rowlen = SCIProwGetNLPNonz(flowrow);
   mcfdata->nnewcols = 0;
   for( i = 0; i < rowlen; i++ )
   {
      SCIP_COL** capacityrowcols;
      int capacityrowlen;
      int arcid;
      int c;
      int j;

      c = SCIPcolGetLPPos(rowcols[i]);
      assert(0 <= c && c < SCIPgetNLPCols(scip));

      /* get arc id of the column in the flow conservation constraint */
      arcid = colarcid[c];
      if( arcid == -1 )
         continue;
      assert(arcid < mcfdata->narcs);

      /* collect flow variables in the capacity constraint of this arc */
      assert(capacityrows[arcid] != NULL);
      capacityrowcols = SCIProwGetCols(capacityrows[arcid]);
      capacityrowlen = SCIProwGetNLPNonz(capacityrows[arcid]);
      for( j = 0; j < capacityrowlen; j++ )
      {
         int caprowc;

         caprowc = SCIPcolGetLPPos(capacityrowcols[j]);
         assert(0 <= caprowc && caprowc < SCIPgetNLPCols(scip));

         /* ignore columns that do not belong to a commodity, i.e., are not flow variables */
         if( colcommodity[caprowc] == -1 )
         {
            assert(colarcid[caprowc] == -1);
            continue;
         }
         assert(colarcid[caprowc] == arcid);

         /* ignore columns in the same commodity as the base row */
         if( colcommodity[caprowc] == basecommodity )
            continue;

         /* if not already done, collect the column */
         if( !colisincident[caprowc] )
         {
            assert(mcfdata->nnewcols < SCIPgetNLPCols(scip));
            colisincident[caprowc] = TRUE;
            newcols[mcfdata->nnewcols] = caprowc;
            mcfdata->nnewcols++;
         }
      }
   }
}

/** compares given row against a base node flow row and calculates a similarity score;
 *  score is 0.0 if the rows are incompatible
 */
static
void getNodeSilimarityScore(
   SCIP*                 scip,               /**< SCIP data structure */ 
   MCFDATA*              mcfdata,            /**< internal MCF extraction data to pass to subroutines */
   int                   baserowlen,         /**< length of base node flow row */
   int*                  basearcpattern,     /**< arc patern of base node flow row */
   SCIP_ROW*             row,                /**< row to compare against base node flow row */
   SCIP_Real*            score,              /**< pointer to store the similarity score */
   SCIP_Bool*            invertcommodity     /**< pointer to store whether the arcs in the commodity of the row have
                                              *   to be inverted for the row to be compatible to the base row */
   )
{
   int* commoditysigns = mcfdata->commoditysigns;
   int* rowcommodity   = mcfdata->rowcommodity;
   int* colarcid       = mcfdata->colarcid;

   SCIP_COL** rowcols;
   SCIP_Real* rowvals;
   int rowlen;
   int rowcom;
   int rowcomsign;
   SCIP_Bool incompatable;
   int overlap;
   int r;
   int i;

   *score = 0.0;
   *invertcommodity = FALSE;

   r = SCIProwGetLPPos(row);
   assert(0 <= r && r < SCIPgetNLPRows(scip));
   rowcom = rowcommodity[r];
   assert(0 <= rowcom && rowcom < mcfdata->ncommodities);
   rowcomsign = commoditysigns[rowcom];
   assert(-1 <= rowcomsign && rowcomsign <= +1);

   rowcols = SCIProwGetCols(row);
   rowvals = SCIProwGetVals(row);
   rowlen = SCIProwGetNLPNonz(row);
   incompatable = FALSE;
   overlap = 0;
   for( i = 0; i < rowlen; i++ )
   {
      int c;
      int arcid;
      int valsign;

      c = SCIPcolGetLPPos(rowcols[i]);
      assert(0 <= c && c < SCIPgetNLPCols(scip));

      arcid = colarcid[c];
      if( arcid == -1 )
         continue;
      assert(arcid < mcfdata->narcs);

      valsign = (rowvals[i] > 0.0 ? +1 : -1);
      if( basearcpattern[arcid] == valsign )
      {
         /* both rows have the same sign for this arc */
         if( rowcomsign == -1 )
         {
            /* commodity is inverted, but we have the same sign: this is incompatible */
            incompatable = TRUE;
            break;
         }
         else
         {
            /* we must not invert arc directions in the commodity */
            rowcomsign = +1;
            overlap++;
         }
      }
      else if( basearcpattern[arcid] == -valsign )
      {
         /* both rows have opposite sign for this arc */
         if( rowcomsign == +1 )
         {
            /* commodity cannot be inverted, but we have opposite signs: this is incompatible */
            incompatable = TRUE;
            break;
         }
         else
         {
            /* we must invert arc directions in the commodity */
            rowcomsign = -1;
            overlap++;
         }
      }
      else
         assert(basearcpattern[arcid] == 0);
   }

   /* calculate the score: maximize overlap and use minimal number of non-overlapping entries as tie breaker */
   if( !incompatable && overlap > 0 )
   {
      assert(overlap <= rowlen);
      assert(overlap <= baserowlen);
      *score = overlap - (rowlen + baserowlen - 2*overlap)/(mcfdata->narcs+1.0);
      *invertcommodity = (rowcomsign == -1);
   }
}


/** assigns node ids to flow conservation constraints */
static
SCIP_RETCODE extractNodes(
   SCIP*                 scip,               /**< SCIP data structure */ 
   MCFDATA*              mcfdata             /**< internal MCF extraction data to pass to subroutines */
   )
{
   unsigned char* flowrowsigns   = mcfdata->flowrowsigns;
   int*           commoditysigns = mcfdata->commoditysigns;
   int*           colcommodity   = mcfdata->colcommodity;
   int*           rowcommodity   = mcfdata->rowcommodity;
   int*           colarcid       = mcfdata->colarcid;
   int*           rowarcid       = mcfdata->rowarcid;
   int*           newcols        = mcfdata->newcols;
   int*           rownodeid;
   SCIP_Bool*     colisincident;
   SCIP_Bool*     rowprocessed;

   SCIP_ROW** rows;
   SCIP_COL** cols;
   int nrows;
   int ncols;

   int* arcpattern;
   SCIP_ROW** bestflowrows;
   SCIP_Real* bestscores;
   SCIP_Bool* bestinverted;
   int r;
   int c;

   assert(mcfdata->nnodes == 0);

   /* get LP data */
   SCIP_CALL( SCIPgetLPRowsData(scip, &rows, &nrows) );
   SCIP_CALL( SCIPgetLPColsData(scip, &cols, &ncols) );

   /* allocate temporary memory */
   SCIP_CALL( SCIPallocMemoryArray(scip, &mcfdata->rownodeid, nrows) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &mcfdata->colisincident, ncols) );
   rownodeid = mcfdata->rownodeid;
   colisincident = mcfdata->colisincident;

   /* allocate temporary local memory */
   SCIP_CALL( SCIPallocMemoryArray(scip, &arcpattern, mcfdata->narcs) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &bestflowrows, mcfdata->ncommodities) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &bestscores, mcfdata->ncommodities) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &bestinverted, mcfdata->ncommodities) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &rowprocessed, nrows) );

   /* initialize temporary memory */
   for( r = 0; r < nrows; r++ )
      rownodeid[r] = -1;
   for( c = 0; c < ncols; c++ )
      colisincident[c] = FALSE;

   /* process all flow conservation constraints that have been used */
   for( r = 0; r < nrows; r++ )
   {
      SCIP_COL** rowcols;
      SCIP_Real* rowvals;
      int rowlen;
      int rowscale;
      int basecommodity;
      int i;

      /* ignore rows that are not used as flow conservation constraint */
      basecommodity = rowcommodity[r];
      if( basecommodity == -1 )
         continue;
      assert((flowrowsigns[r] & (LHSASSIGNED | RHSASSIGNED)) != 0);
      assert(rowarcid[r] == -1);

      /* skip rows that are already assigned to a node */
      if( rownodeid[r] >= 0 )
         continue;

      /* assign row to new node id */
      SCIPdebugMessage("assigning row %d <%s> of commodity %d to node %d\n", r, SCIProwGetName(rows[r]), basecommodity, mcfdata->nnodes);
      rownodeid[r] = mcfdata->nnodes;

      /* get the arc pattern of the flow row */
      BMSclearMemoryArray(arcpattern, mcfdata->narcs);
      rowcols = SCIProwGetCols(rows[r]);
      rowvals = SCIProwGetVals(rows[r]);
      rowlen = SCIProwGetNLPNonz(rows[r]);
      if( (flowrowsigns[r] & LHSASSIGNED) != 0 )
         rowscale = +1;
      else
         rowscale = -1;
      for( i = 0; i < rowlen; i++ )
      {
         int arcid;

         c = SCIPcolGetLPPos(rowcols[i]);
         assert(0 <= c && c < ncols);
         arcid = colarcid[c];
         if( arcid >= 0 )
         {
            if( rowvals[i] > 0.0 )
               arcpattern[arcid] = rowscale;
            else
               arcpattern[arcid] = -rowscale;
         }
      }

      /* initialize arrays to store best flow rows */
      for( i = 0; i < mcfdata->ncommodities; i++ )
      {
         bestflowrows[i] = NULL;
         bestscores[i] = 0.0;
         bestinverted[i] = FALSE;
      }

      /* collect columns that are member of incident arc capacity constraints */
      collectIncidentFlowCols(scip, mcfdata, rows[r], basecommodity);

      /* initialize rowprocessed array */
      BMSclearMemoryArray(rowprocessed, nrows);

      /* identify flow conservation constraints in other commodities that match this node;
       * search for flow rows in the column vectors of the indicent columns
       */
      for( i = 0; i < mcfdata->nnewcols; i++ )
      {
         SCIP_ROW** colrows;
         int collen;
         int j;

         c = newcols[i];
         assert(0 <= c && c < ncols);
         assert(colcommodity[c] >= 0);
         assert(colcommodity[c] != basecommodity);

         /* clean up the marker array */
         assert(colisincident[c]);
         colisincident[c] = FALSE;

         /* scan column vector for flow conservation constraints */
         colrows = SCIPcolGetRows(cols[c]);
         collen = SCIPcolGetNLPNonz(cols[c]);
         for( j = 0; j < collen; j++ )
         {
            int colr;
            int rowcom;
            SCIP_Real score;
            SCIP_Bool invertcommodity;

            colr = SCIProwGetLPPos(colrows[j]);
            assert(0 <= colr && colr < nrows);

            /* ignore rows that have already been processed */
            if( rowprocessed[colr] )
               continue;
            rowprocessed[colr] = TRUE;

            /* ignore rows that are not flow conservation constraints in the network */
            rowcom = rowcommodity[colr];
            assert(rowcom != basecommodity);
            if( rowcom == -1 )
               continue;
            assert(rowcom == colcommodity[c]);
            assert((flowrowsigns[colr] & (LHSASSIGNED | RHSASSIGNED)) != 0);
            assert(rowarcid[colr] == -1);

            /* ignore rows that are already assigned to a node */
            if( rownodeid[colr] >= 0 )
               continue;

            /* compare row against arc pattern and calculate score */
            getNodeSilimarityScore(scip, mcfdata, rowlen, arcpattern, colrows[j], &score, &invertcommodity);
            if( score > bestscores[rowcom] )
            {
               bestflowrows[rowcom] = colrows[j];
               bestscores[rowcom] = score;
               bestinverted[rowcom] = invertcommodity;
            }
         }
      }

      /* for each commodity, pick the best flow conservation constraint to define this node */
      for( i = 0; i < mcfdata->ncommodities; i++ )
      {
         int comr;

         if( bestflowrows[i] == NULL )
            continue;

         comr = SCIProwGetLPPos(bestflowrows[i]);
         assert(0 <= comr && comr < nrows);
         assert(rowcommodity[comr] == i);
         assert((flowrowsigns[comr] & (LHSASSIGNED | RHSASSIGNED)) != 0);
         assert(rownodeid[comr] == -1);

         /* assign flow row to current node */
         SCIPdebugMessage(" -> assigning row %d <%s> of commodity %d to node %d\n", comr, SCIProwGetName(rows[comr]), i, mcfdata->nnodes);
         rownodeid[comr] = mcfdata->nnodes;

         /* fix the direction of the arcs of the commodity */
         if( bestinverted[i] )
         {
            assert(commoditysigns[i] != +1);
            commoditysigns[i] = -1;
         }
         else
         {
            assert(commoditysigns[i] != -1);
            commoditysigns[i] = +1;
         }
      }
      
      /* increase number of nodes */
      mcfdata->nnodes++;
   }

   /* free local temporary memory */
   SCIPfreeMemoryArray(scip, &rowprocessed);
   SCIPfreeMemoryArray(scip, &bestinverted);
   SCIPfreeMemoryArray(scip, &bestscores);
   SCIPfreeMemoryArray(scip, &bestflowrows);
   SCIPfreeMemoryArray(scip, &arcpattern);

   return SCIP_OKAY;
}

/* if there are still undecided commodity signs, fix them to +1 */
static
void fixCommoditySigns(
   SCIP*                 scip,               /**< SCIP data structure */ 
   MCFDATA*              mcfdata             /**< internal MCF extraction data to pass to subroutines */
   )
{
   int* commoditysigns = mcfdata->commoditysigns;
   int k;

   for( k = 0; k < mcfdata->ncommodities; k++ )
   {
      if( commoditysigns[k] == 0 )
         commoditysigns[k] = +1;
   }
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
   MCFDATA mcfdata;

   SCIP_ROW** rows;
   SCIP_COL** cols;
   int nrows;
   int ncols;

   assert(mcfnetwork != NULL);

   /* Algorithm to identify multi-commodity-flow network with capacity constraints
    *
    * 1. Identify candidate rows for flow conservation constraints and capacity constraints in the LP.
    * 2. Sort flow conservation and capacity constraint candidates by a ranking on
    *    how sure we are that it is indeed a constraint of the desired type.
    * 3. Extract network structure of flow conservation constraints:
    *    (a) Initialize plusflow[c] = minusflow[c] = FALSE for all columns c and other local data.
    *    (b) As long as there are flow conservation candidates left:
    *        (i) Create new commodity and use first flow conservation constraint as newrow.
    *       (ii) Add newrow to commodity, update pluscom/minuscom accordingly.
    *      (iii) For the newly added columns search for an incident flow conservation constraint. Pick the one of highest ranking.
    *       (iv) If found, set newrow to this row and goto (ii).
    * 4. Identify capacity constraints for the arcs and assign arc ids to columns and capacity constraints.
    * 5. Assign node ids to flow conservation constraints.
    */

   /* create network data structure */
   SCIP_CALL( mcfnetworkCreate(scip, mcfnetwork) );
   assert(*mcfnetwork != NULL);
   mcf = *mcfnetwork;

   /* get LP data */
   SCIP_CALL( SCIPgetLPRowsData(scip, &rows, &nrows) );
   SCIP_CALL( SCIPgetLPColsData(scip, &cols, &ncols) );

   /* initialize local extraction data */
   mcfdata.flowrowsigns = NULL;
   mcfdata.flowrowscalars = NULL;
   mcfdata.flowrowscores = NULL;
   mcfdata.capacityrowsigns = NULL;
   mcfdata.capacityrowscores = NULL;
   mcfdata.flowcands = NULL;
   mcfdata.nflowcands = 0;
   mcfdata.capacitycands = NULL;
   mcfdata.ncapacitycands = 0;
   mcfdata.plusflow = NULL;
   mcfdata.minusflow = NULL;
   mcfdata.ncommodities = 0;
   mcfdata.commoditysigns = NULL;
   mcfdata.commoditysignssize = 0;
   mcfdata.colcommodity = NULL;
   mcfdata.rowcommodity = NULL;
   mcfdata.colarcid = NULL;
   mcfdata.rowarcid = NULL;
   mcfdata.rownodeid = NULL;
   mcfdata.newcols = NULL;
   mcfdata.nnewcols = 0;
   mcfdata.narcs = 0;
   mcfdata.nnodes = 0;
   mcfdata.capacityrows = NULL;
   mcfdata.capacityrowssize = 0;
   mcfdata.colisincident = NULL;

   /* 1. Identify candidate rows for flow conservation constraints and capacity constraints in the LP.
    * 2. Sort flow conservation and capacity constraint candidates by a ranking on
    *    how sure we are that it is indeed a constraint of the desired type.
    */
   SCIP_CALL( extractRows(scip, &mcfdata) );

   if( mcfdata.nflowcands > 0 && mcfdata.ncapacitycands > 0 )
   {
      /* 3. Extract network structure of flow conservation constraints. */
      SCIPdebugMessage("****** extracting flow ******\n");
      SCIP_CALL( extractFlow(scip, &mcfdata) );
#ifdef SCIP_DEBUG
      printCommodities(scip, &mcfdata);
#endif
      
      /* 4. Identify capacity constraints for the arcs and assign arc ids to columns and capacity constraints. */
      SCIPdebugMessage("****** extracting capacities ******\n");
      SCIP_CALL( extractCapacities(scip, &mcfdata) );

      /* 5. Assign node ids to flow conservation constraints. */
      SCIPdebugMessage("****** extracting nodes ******\n");
      SCIP_CALL( extractNodes(scip, &mcfdata) );

      /* if there are still undecided commodity signs, fix them to +1 */
      fixCommoditySigns(scip, &mcfdata);
   }

#ifdef SCIP_DEBUG
   printCommodities(scip, &mcfdata);
#endif

   /* free memory */
   SCIPfreeMemoryArrayNull(scip, &mcfdata.colisincident);
   SCIPfreeMemoryArrayNull(scip, &mcfdata.capacityrows);
   SCIPfreeMemoryArrayNull(scip, &mcfdata.rownodeid);
   SCIPfreeMemoryArrayNull(scip, &mcfdata.rowarcid);
   SCIPfreeMemoryArrayNull(scip, &mcfdata.colarcid);
   SCIPfreeMemoryArrayNull(scip, &mcfdata.newcols);
   SCIPfreeMemoryArrayNull(scip, &mcfdata.rowcommodity);
   SCIPfreeMemoryArrayNull(scip, &mcfdata.colcommodity);
   SCIPfreeMemoryArrayNull(scip, &mcfdata.commoditysigns);
   SCIPfreeMemoryArrayNull(scip, &mcfdata.minusflow);
   SCIPfreeMemoryArrayNull(scip, &mcfdata.plusflow);
   SCIPfreeMemoryArrayNull(scip, &mcfdata.capacitycands);
   SCIPfreeMemoryArrayNull(scip, &mcfdata.flowcands);
   SCIPfreeMemoryArrayNull(scip, &mcfdata.capacityrowscores);
   SCIPfreeMemoryArrayNull(scip, &mcfdata.capacityrowsigns);
   SCIPfreeMemoryArrayNull(scip, &mcfdata.flowrowscores);
   SCIPfreeMemoryArrayNull(scip, &mcfdata.flowrowscalars);
   SCIPfreeMemoryArrayNull(scip, &mcfdata.flowrowsigns);

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
