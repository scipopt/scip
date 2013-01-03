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

/**@file   sepa_flowcover.c
 * @brief  flow cover cuts separator
 * @author Kati Wolter
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>

#include "scip/sepa_flowcover.h"
#include "scip/cons_knapsack.h"
#include "scip/pub_misc.h"


#define SEPA_NAME              "flowcover"
#define SEPA_DESC              "flow cover cuts separator (c-MIR approach)"
#define SEPA_PRIORITY             -4000
#define SEPA_FREQ                     0
#define SEPA_MAXBOUNDDIST           0.0
#define SEPA_USESSUBSCIP          FALSE /**< does the separator use a secondary SCIP instance? */
#define SEPA_DELAY                FALSE /**< should separation method be delayed, if other separators found cuts? */

#define DEFAULT_MAXROUNDS             5 /**< maximal number of separation rounds per node (-1: unlimited) */
#define DEFAULT_MAXROUNDSROOT        15 /**< maximal number of separation rounds in the root node (-1: unlimited) */
#define DEFAULT_MAXTRIES            100 /**< maximal number of rows to separate flow cover cuts for per separation round 
                                         *   (-1: unlimited) */
#define DEFAULT_MAXTRIESROOT         -1 /**< maximal number of rows to separate flow cover cuts for per separation round 
                                         *   in the root (-1: unlimited) */
#define DEFAULT_MAXFAILS             50 /**< maximal number of consecutive fails to generate a cut per separation round
                                        *   (-1: unlimited) */
#define DEFAULT_MAXFAILSROOT        100 /**< maximal number of consecutive fails to generate a cut per separation round
                                         *   in the root (-1: unlimited) */
#define DEFAULT_MAXSEPACUTS         100 /**< maximal number of flow cover cuts separated per separation round */
#define DEFAULT_MAXSEPACUTSROOT     200 /**< maximal number of flow cover cuts separated per separation round in the root */
#define DEFAULT_MAXSLACK  SCIP_REAL_MAX /**< maximal slack of rows to separate flow cover cuts for */
#define DEFAULT_MAXSLACKROOT SCIP_REAL_MAX /**< maximal slack of rows to separate flow cover cuts for in the root */
#define DEFAULT_SLACKSCORE        1e-03 /**< weight of slack in the scoring of the rows */
#define DEFAULT_MAXROWDENSITY       1.0 /**< maximal density of rows to separate flow cover cuts for */
#define DEFAULT_DYNAMICCUTS        TRUE /**< should generated cuts be removed from the LP if they are no longer tight? */
#define DEFAULT_MAXTESTDELTA         10 /**< cut generation heuristic: maximal number of different deltas to try */
#define DEFAULT_MULTBYMINUSONE     TRUE /**< should flow cover cuts be separated for 0-1 single node flow set with reversed arcs in addition? */

#define BOUNDSWITCH                 0.5
#define ALLOWLOCAL                 TRUE 
#define DENSSCORE                 1e-04
/*#define MAKECONTINTEGRAL          FALSE*/
#define MINFRAC                    0.01
#define MAXFRAC                    0.95
#define FIXINTEGRALRHS            FALSE
#define MAXDNOM                  1000LL 
#define MINDELTA                  1e-03 
#define MAXDELTA                  1e-09 
#define MAXSCALE                 1000.0 
#define MAXDYNPROGSPACE         1000000 

#define MAXAGGRLEN(nvars)          (0.1*(nvars)+1000) /**< maximal length of base inequality */
#define MAXABSVBCOEF               1e+5 /**< maximal absolute coefficient in variable bounds used for snf relaxation */
#define MAXBOUND                  1e+10 /**< maximal value of normal bounds used for snf relaxation */


/*
 * Data structures
 */

/** separator data */
struct SCIP_SepaData
{
   int                   maxrounds;          /**< maximal number of separation rounds per node (-1: unlimited) */
   int                   maxroundsroot;      /**< maximal number of separation rounds in the root node (-1: unlimited) */
   int                   maxtries;           /**< maximal number of rows to separate flow cover cuts for per separation round 
                                              *   (-1: unlimited) */
   int                   maxtriesroot;       /**< maximal number of rows to separate flow cover cuts for per separation round 
                                              *   in the root (-1: unlimited) */
   int                   maxfails;           /**< maximal number of consecutive fails to generate a cut per separation round
                                              *   (-1: unlimited) */
   int                   maxfailsroot;       /**< maximal number of consecutive fails to generate a cut per separation round
                                              *   in the root (-1: unlimited) */
   int                   maxsepacuts;        /**< maximal number of flow cover cuts separated per separation round */
   int                   maxsepacutsroot;    /**< maximal number of flow cover cuts separated per separation round in the root */
   SCIP_Real             maxslack;           /**< maximal slack of rows to separate flow cover cuts for */
   SCIP_Real             maxslackroot;       /**< maximal slack of rows to separate flow cover cuts for in the root */
   SCIP_Real             slackscore;         /**< weight of slack in the scoring of the rows */
   SCIP_Real             maxrowdensity;      /**< maximal density of rows to separate flow cover cuts for */
   SCIP_Bool             dynamiccuts;        /**< should generated cuts be removed from the LP if they are no longer tight? */
   SCIP_Bool             multbyminusone;     /**< should flow cover cuts be separated for 0-1 single node flow set with reversed arcs in addition? */
   int                   maxtestdelta;       /**< cut generation heuristic: maximal number of different deltas to try */
};


/*
 * Local methods
 */

/** get LP solution value and index of variable lower bound (with binary variable) which is closest to the current LP 
 *  solution value of a given variable; candidates have to meet certain criteria in order to ensure the nonnegativity 
 *  of the variable upper bound imposed on the real variable in the 0-1 single node flow relaxation associated with the
 *  given variable
 */
static
SCIP_RETCODE getClosestVlb(
   SCIP*                 scip,               /**< SCIP data structure */ 
   SCIP_VAR*             var,                /**< given active problem variable */
   SCIP_Real             bestsub,            /**< closest simple upper bound of given variable */
   SCIP_Real             rowcoef,            /**< coefficient of given variable in current row */
   SCIP_Real*            rowcoefsbinary,     /**< coefficient of all binary problem variables in current row */ 
   SCIP_Real*            varsolvals,         /**< LP solution value of all active problem variables */
   int*                  assoctransvars,     /**< associated var in relaxed set for all vars of row; construction is not finished yet */ 
   SCIP_Real*            closestvlb,         /**< pointer to store the LP sol value of the closest variable lower bound */
   int*                  closestvlbidx       /**< pointer to store the index of the closest vlb; -1 if no vlb was found */
   )
{
   int nvlbs;

   assert(scip != NULL);
   assert(var != NULL);
   assert(bestsub == SCIPvarGetUbGlobal(var) || bestsub == SCIPvarGetUbLocal(var)); /*lint !e777*/
   assert(!SCIPisInfinity(scip, bestsub));
   assert(!SCIPisZero(scip, rowcoef));
   assert(rowcoefsbinary != NULL);
   assert(varsolvals != NULL);
   assert(assoctransvars != NULL);
   assert(closestvlb != NULL);
   assert(closestvlbidx != NULL);

   nvlbs = SCIPvarGetNVlbs(var);

   *closestvlbidx = -1;
   *closestvlb = -SCIPinfinity(scip); 
   if( nvlbs > 0 )
   {
      SCIP_VAR** vlbvars;
      SCIP_Real* vlbcoefs;
      SCIP_Real* vlbconsts;
      int i;
      
      vlbvars = SCIPvarGetVlbVars(var);
      vlbcoefs = SCIPvarGetVlbCoefs(var);
      vlbconsts = SCIPvarGetVlbConstants(var);

      for( i = 0; i < nvlbs; i++ )
      {
         SCIP_Real rowcoefbinvar;
         SCIP_Real val1;
         SCIP_Real val2;
         SCIP_Real vlbsol;
         SCIP_Bool meetscriteria; 
         int probidxbinvar;
      
         /* use only variable lower bounds l~_i * x_i + d_i with x_i binary which are active */
         if( !SCIPvarIsBinary(vlbvars[i])  || !SCIPvarIsActive(vlbvars[i]) )
            continue;
         
         /* check if current variable lower bound l~_i * x_i + d_i imposed on y_j meets the following criteria:
          * (let a_j  = coefficient of y_j in current row,
          *      u_j  = closest simple upper bound imposed on y_j,
          *      c_i  = coefficient of x_i in current row)
          *   0. no other non-binary variable y_k has used a variable bound with x_i to get transformed variable y'_k yet
          * if a_j > 0: 
          *   1. u_j <= d_i
          *   2. a_j ( u_j - d_i ) + c_i <= 0
          *   3. a_j l~_i + c_i <= 0
          * if a_j < 0: 
          *   1. u_j <= d_i
          *   2. a_j ( u_j - d_i ) + c_i >= 0
          *   3. a_j l~_i + c_i >= 0 
          */
         probidxbinvar = SCIPvarGetProbindex(vlbvars[i]);
         rowcoefbinvar = rowcoefsbinary[probidxbinvar];

         val1 = ( rowcoef * ( bestsub - vlbconsts[i] ) ) + rowcoefbinvar;
         val2 = ( rowcoef * vlbcoefs[i] ) + rowcoefbinvar;

         meetscriteria = FALSE;
         if( SCIPisPositive(scip, rowcoef) )
         {
            if( assoctransvars[probidxbinvar] == -1 && SCIPisFeasLE(scip, bestsub, vlbconsts[i])
               && SCIPisFeasLE(scip, val1, 0.0) && SCIPisFeasLE(scip, val2, 0.0) )
               meetscriteria = TRUE;
         }
         else
         {
            assert(SCIPisNegative(scip, rowcoef));
            if( assoctransvars[probidxbinvar] == -1 && SCIPisFeasLE(scip, bestsub, vlbconsts[i])
               && SCIPisFeasGE(scip, val1, 0.0) && SCIPisFeasGE(scip, val2, 0.0) )
               meetscriteria = TRUE;
         }
   
         /* variable lower bound does not meet criteria */
         if( !meetscriteria )
            continue;

         /* for numerical reasons, ignore variable bounds with large absolute coefficient and 
          * those which lead to an infinite variable bound coefficient (val2) in snf relaxation 
          */
         if( REALABS(vlbcoefs[i]) > MAXABSVBCOEF || SCIPisInfinity(scip, REALABS(val2)) )
            continue;

         vlbsol = (vlbcoefs[i] * varsolvals[probidxbinvar]) + vlbconsts[i];
         if( SCIPisGT(scip, vlbsol, *closestvlb) )
         {
            *closestvlb = vlbsol;
            *closestvlbidx = i;
         }
         assert(*closestvlbidx >= 0);

      }
   }

   return SCIP_OKAY;
}

/** get LP solution value and index of variable upper bound (with binary variable) which is closest to the current LP 
 *  solution value of a given variable; candidates have to meet certain criteria in order to ensure the nonnegativity 
 *  of the variable upper bound imposed on the real variable in the 0-1 single node flow relaxation associated with the
 *  given variable
 */
static
SCIP_RETCODE getClosestVub(
   SCIP*                 scip,               /**< SCIP data structure */ 
   SCIP_VAR*             var,                /**< given active problem variable */
   SCIP_Real             bestslb,            /**< closest simple lower bound of given variable */
   SCIP_Real             rowcoef,            /**< coefficient of given variable in current row */
   SCIP_Real*            rowcoefsbinary,     /**< coefficient of all binary problem variables in current row */ 
   SCIP_Real*            varsolvals,         /**< LP solution value of all active problem variables */
   int*                  assoctransvars,     /**< associated var in relaxed set for all vars of row; construction is not finished yet */ 
   SCIP_Real*            closestvub,         /**< pointer to store the LP sol value of the closest variable upper bound */
   int*                  closestvubidx       /**< pointer to store the index of the closest vub; -1 if no vub was found */
   )
{
   int nvubs;

   assert(scip != NULL);
   assert(var != NULL);
   assert(bestslb == SCIPvarGetLbGlobal(var) || bestslb == SCIPvarGetLbLocal(var)); /*lint !e777*/
   assert(!SCIPisInfinity(scip, - bestslb));
   assert(!SCIPisZero(scip, rowcoef));
   assert(rowcoefsbinary != NULL);
   assert(varsolvals != NULL);
   assert(assoctransvars != NULL);
   assert(closestvub != NULL);
   assert(closestvubidx != NULL);

   nvubs = SCIPvarGetNVubs(var);

   *closestvubidx = -1;
   *closestvub = SCIPinfinity(scip);
   if( nvubs > 0 )
   {
      SCIP_VAR** vubvars;
      SCIP_Real* vubcoefs;
      SCIP_Real* vubconsts;
      int i;

      vubvars = SCIPvarGetVubVars(var);
      vubcoefs = SCIPvarGetVubCoefs(var);
      vubconsts = SCIPvarGetVubConstants(var);
      
      for( i = 0; i < nvubs; i++ )
      {
         SCIP_Real rowcoefbinvar;
         SCIP_Real val1;
         SCIP_Real val2;
         SCIP_Real vubsol;
         SCIP_Bool meetscriteria; 
         int probidxbinvar;
    
         /* use only variable upper bound u~_i * x_i + d_i with x_i binary and which are active */
         if( !SCIPvarIsBinary(vubvars[i]) || !SCIPvarIsActive(vubvars[i]))
            continue;

         /* checks if current variable upper bound u~_i * x_i + d_i meets the following criteria
          * (let a_j  = coefficient of y_j in current row,
          *      l_j  = closest simple lower bound imposed on y_j,
          *      c_i  = coefficient of x_i in current row)
          *   0. no other non-binary variable y_k has used a variable bound with x_i to get transformed variable y'_k 
          * if a > 0: 
          *   1. l_j >= d_i
          *   2. a_j ( l_i - d_i ) + c_i >= 0
          *   3. a_j u~_i + c_i >= 0
          * if a < 0: 
          *   1. l_j >= d_i
          *   2. a_j ( l_j - d_i ) + c_i <= 0
          *   3. a_j u~_i + c_i <= 0 
          */
         probidxbinvar = SCIPvarGetProbindex(vubvars[i]);
         rowcoefbinvar = rowcoefsbinary[probidxbinvar];

         val1 = ( rowcoef * ( bestslb - vubconsts[i] ) ) + rowcoefbinvar;
         val2 = ( rowcoef * vubcoefs[i] ) + rowcoefbinvar;

         meetscriteria = FALSE;
         if( SCIPisPositive(scip, rowcoef) )
         {
            if( assoctransvars[probidxbinvar] == -1 && SCIPisFeasGE(scip, bestslb, vubconsts[i])
               && SCIPisFeasGE(scip, val1, 0.0) && SCIPisFeasGE(scip, val2, 0.0) )
               meetscriteria = TRUE;
         }
         else
         {
            assert(SCIPisNegative(scip, rowcoef));
            if( assoctransvars[probidxbinvar] == -1 && SCIPisFeasGE(scip, bestslb, vubconsts[i])
               && SCIPisFeasLE(scip, val1, 0.0) && SCIPisFeasLE(scip, val2, 0.0) )
               meetscriteria = TRUE;
         }

         /* variable upper bound does not meet criteria */
         if( !meetscriteria )
            continue;
        
         /* for numerical reasons, ignore variable bounds with large absolute coefficient and
          * those which lead to an infinite variable bound coefficient (val2) in snf relaxation 
          */
         if( REALABS(vubcoefs[i]) > MAXABSVBCOEF || SCIPisInfinity(scip, REALABS(val2)) )
            continue;

         vubsol = vubcoefs[i] * varsolvals[probidxbinvar] + vubconsts[i];
         if( SCIPisLT(scip, vubsol, *closestvub) )
         {
            *closestvub = vubsol;
            *closestvubidx = i;
         }
         assert(*closestvubidx >= 0);
      }
   }

   return SCIP_OKAY;
}

/** return global or local lower bound of given variable whichever is closer to the variables current LP solution value */
static
void getClosestLb(
   SCIP*                 scip,               /**< SCIP data structure */ 
   SCIP_VAR*             var,                /**< active problem variable */
   SCIP_Bool             allowlocal,         /**< should local information allowed to be used, resulting in a local cut? */
   SCIP_Real*            closestlb,          /**< pointer to store the value of the closest variable lower bound */
   int*                  closestlbtype       /**< pointer to store type of closest bound; -1 if global lb, -2 otherwise */
   )
{
   assert(closestlb != NULL);
   assert(closestlbtype != NULL);
   
   *closestlb = SCIPvarGetLbGlobal(var);
   *closestlbtype = -1;
   if( allowlocal )
   {
      SCIP_Real loclb;
      loclb = SCIPvarGetLbLocal(var);
      if( SCIPisGT(scip, loclb, *closestlb) )
      {
         *closestlb = loclb;
         *closestlbtype = -2;
      }
   }

   /* due to numerical reasons, huge bounds are relaxed to infinite bounds; this way the bounds are not used for
    * the construction of the 0-1 single node flow relaxation
    */
   if( *closestlb <= -MAXBOUND )
      *closestlb = -SCIPinfinity(scip);
}

/** return global or local upper bound of given variable whichever is closer to the variables current LP solution value */
static
void getClosestUb(
   SCIP*                 scip,               /**< SCIP data structure */ 
   SCIP_VAR*             var,                /**< active problem variable */
   SCIP_Bool             allowlocal,         /**< should local information allowed to be used, resulting in a local cut? */
   SCIP_Real*            closestub,          /**< pointer to store the value of the closest upper bound */
   int*                  closestubtype       /**< pointer to store type of closest bound; -1 if global ub, -2 otherwise */
   )
{
   assert(closestub != NULL);
   assert(closestubtype != NULL);
   
   *closestub = SCIPvarGetUbGlobal(var);
   *closestubtype = -1;
   if( allowlocal )
   {
      SCIP_Real locub;
      locub = SCIPvarGetUbLocal(var);
      if( SCIPisLT(scip, locub, *closestub) )
      {
         *closestub = locub;
         *closestubtype = -2;
      }
   }

   /* due to numerical reasons, huge bounds are relaxed to infinite bounds; this way the bounds are not used for
    * the construction of the 0-1 single node flow relaxation
    */
   if( *closestub >= MAXBOUND )
      *closestub = SCIPinfinity(scip);
}

/** construct a 0-1 single node flow relaxation (with some additional simple constraints) of a mixed integer set 
 *  corresponding to the given row lhs <= a * x + const <= rhs; depending on the given values rowweight and scale
 *  the mixed integer set which should be used is defined by the mixed integer constraint  
 *    a * (x,y)     <=    rhs - const      if (rowweight =  1, scale =  1) or (rowweight = -1, scale = -1, rhs <  infinity)
 *  - a * (x,y)     <= - (lhs - const)     if (rowweight = -1, scale =  1) or (rowweight =  1, scale = -1, lhs > -infinity)
 *  - a * (x,y) - s <= - (rhs - const)     if (rowweight =  1, scale = -1, lhs = -infinity)    
 *    a * (x,y) - s <=   (lhs - const)     if (rowweight = -1, scale = -1, rhs =  infinity)
 */
static
SCIP_RETCODE constructSNFRelaxation(
   SCIP*                 scip,               /**< SCIP data structure */ 
   SCIP_VAR**            vars,               /**< active problem variables */
   int                   nvars,              /**< number of active problem variables */
   SCIP_Real*            varsolvals,         /**< solution values of active problem variables */
   SCIP_ROW*             row,                /**< given row */
   SCIP_Real             rowweight,          /**< weight of given row; can be +1 or -1 */
   SCIP_Real             scale,              /**< additional scaling factor for given row */
   int*                  boundsfortrans,     /**< pointer to store bound used for all non-binary vars of row */
   SCIP_BOUNDTYPE*       boundtypesfortrans, /**< pointer to store type of bound used for all non-binary vars of row */
   int*                  assoctransvars,     /**< pointer to store associated var in relaxed set for all vars of row */ 
   int*                  transvarcoefs,      /**< pointer to store coefficient of all vars in relaxed set */ 
   SCIP_Real*            transbinvarsolvals, /**< pointer to store sol val of bin var in vub of all vars in relaxed set */
   SCIP_Real*            transcontvarsolvals,/**< pointer to store sol val of all real vars in relaxed set */
   SCIP_Real*            transvarvubcoefs,   /**< pointer to store coefficient in vub of all vars in relaxed set */
   int*                  ntransvars,         /**< pointer to store number of vars in relaxed set */
   SCIP_Real*            transrhs,           /**< pointer to store rhs in relaxed set */ 
   SCIP_Bool*            success             /**< pointer to store whether a relaxation was constructed */
   )
{
   SCIP_COL** nonzcols;   
   SCIP_Real* nonzcoefs;  
   SCIP_Real* rowcoefsbinary;  
   int* nonzcolsbinary;
   int* nonzcolsnonbinary;
   int nnonzcols;         
   int nnonzcolsbinary;
   int nnonzcolsnonbinary;
   int c;

   assert(scip != NULL);
   assert(vars!= NULL);
   assert(nvars > 0);
   assert(varsolvals != NULL);
   assert(row != NULL);
   assert(rowweight == 1.0 || rowweight == -1.0); 
   assert(( rowweight == 1.0 && !SCIPisInfinity(scip, SCIProwGetRhs(row)) ) 
      || ( rowweight == -1.0 && !SCIPisInfinity(scip, -SCIProwGetLhs(row)) ));
   assert(scale == 1.0 || scale == -1.0);
   assert(boundsfortrans != NULL);
   assert(boundtypesfortrans != NULL);
   assert(assoctransvars != NULL);
   assert(transvarcoefs != NULL);
   assert(transbinvarsolvals != NULL);
   assert(transcontvarsolvals != NULL);
   assert(transvarvubcoefs != NULL);
   assert(ntransvars != NULL);
   assert(transrhs != NULL);
   assert(success != NULL);

   *success = FALSE;

   SCIPdebugMessage("--------------------- construction of SNF relaxation ------------------------------------\n");
   
   /* get nonzero columns and coefficients of row */
   nonzcols =  SCIProwGetCols(row);
   nnonzcols = SCIProwGetNLPNonz(row);
   nonzcoefs = SCIProwGetVals(row);
   SCIPdebugMessage("nnonzcols = %d\n",nnonzcols);

   /* get data structures */
   SCIP_CALL( SCIPallocBufferArray(scip, &nonzcolsbinary, nnonzcols) );
   SCIP_CALL( SCIPallocBufferArray(scip, &nonzcolsnonbinary, nnonzcols) );
   SCIP_CALL( SCIPallocBufferArray(scip, &rowcoefsbinary, nvars) );

   /* store nonzero columns representing binary and non-binary variables, and get active binary problem variables the
    * coefficient in the row 
    */
   nnonzcolsbinary = 0;
   nnonzcolsnonbinary = 0;
   BMSclearMemoryArray(rowcoefsbinary, nvars);
   for( c = 0; c < nnonzcols; c++ )
   {
      SCIP_COL* col; 
      SCIP_VAR* var;

      col = nonzcols[c];
      var = SCIPcolGetVar(col);

      assert(!SCIPisZero(scip, nonzcoefs[c]));

      if( SCIPvarIsBinary(var) )
      {
         /* saves column for binary variable */
         nonzcolsbinary[nnonzcolsbinary] = c;
         nnonzcolsbinary++;

         /* saves row coefficient of binary variable */ 
         assert(SCIPvarGetProbindex(var) > -1 &&  SCIPvarGetProbindex(var) < nvars);
         rowcoefsbinary[SCIPvarGetProbindex(var)] = rowweight * scale * nonzcoefs[c];
      }
      else
      {
         /* saves column for non-binary variable */
         nonzcolsnonbinary[nnonzcolsnonbinary] = c;
         nnonzcolsnonbinary++;
      }
   }
   assert(nnonzcolsbinary + nnonzcolsnonbinary == nnonzcols); 

   /* initialize data structures */
   for( c = 0; c < nvars; c++ )
   {
      assoctransvars[c] = -1;
      boundsfortrans[c] = -3;
   }
   *ntransvars = 0;

   /* initialize right hand side of constraint in 0-1 single node flow relaxation */
   if( rowweight * scale == 1.0 && !SCIPisInfinity(scip, SCIProwGetRhs(row)) )
      *transrhs = SCIProwGetRhs(row) - SCIProwGetConstant(row);
   else if( rowweight * scale == 1.0 && SCIPisInfinity(scip, SCIProwGetRhs(row)) )
   {
      assert(rowweight == -1.0 && scale == -1.0);
      *transrhs = SCIProwGetLhs(row) - SCIProwGetConstant(row);
   }
   else if( rowweight * scale == -1.0 && !SCIPisInfinity(scip, -SCIProwGetLhs(row)) )
      *transrhs = - SCIProwGetLhs(row) + SCIProwGetConstant(row); 
   else
   {
      assert(rowweight == 1.0 && scale == -1.0 && SCIPisInfinity(scip, -SCIProwGetLhs(row)));
      *transrhs = - SCIProwGetRhs(row) + SCIProwGetConstant(row);
   }
   
   /* for each non-binary variable y_j in the row with nonzero row coefficient perform
    *   1. get closest simple or variable lower bound and closest simple or variable upper bound
    *   2. decide which bound is used to define the real variable y'_j in the 0-1 single node flow relaxation
    *   3. construct y'_j with 0 <= y'_j <= u'_j x_j
    *   4. store for y_j and x_j (if x_j is a binary variable in the row) that y'_j is the associated real variable 
    *      in the 0-1 single node flow relaxation and for y_j the bound used to define y'_j.
    *
    * for each binary variable x_j in the row which has not been handled with a non-binary variable perform
    *   1. construct y'_j with 0 <= y'_j <= u'_j x_j
    *   2. store for x_j that y'_j is the associated real variable in the 0-1 single node flow relaxation. 
    *  
    * start with non-binary variables because a binary variable x_j which is involved in a used variable bound 
    * imposed on a non-binary variable y_j has to be handled together with the non-binary variable y_j. 
    */
   SCIPdebugMessage("transformation for NONBINARY variables (nnonbinvars=%d):\n", nnonzcolsnonbinary);

   /* non-binary variables and binary variables contained in used variable bounds */
   for( c = 0; c < nnonzcolsnonbinary; c++ )
   {
      SCIP_VAR* var;
      SCIP_Real bestlb;
      SCIP_Real bestub;
      SCIP_Real bestslb;
      SCIP_Real bestsub;
      SCIP_Real rowcoef; 
      SCIP_Bool uselb;
      int bestlbtype;
      int bestubtype;
      int bestslbtype;
      int bestsubtype;
      int probidx;
      
      bestlb = -SCIPinfinity(scip);
      bestub = SCIPinfinity(scip);
      bestlbtype = -3;
      bestubtype = -3;

      var = SCIPcolGetVar(nonzcols[nonzcolsnonbinary[c]]);
      probidx = SCIPvarGetProbindex(var);
      rowcoef = rowweight * scale * nonzcoefs[nonzcolsnonbinary[c]];

      assert(assoctransvars[probidx] == -1);
      assert(boundsfortrans[probidx] == -3);
      assert(!SCIPisZero(scip, rowcoef));

      /* get closest simple lower bound and closest simple upper bound */
      getClosestLb(scip, var, ALLOWLOCAL, &bestslb, &bestslbtype);
      getClosestUb(scip, var, ALLOWLOCAL, &bestsub, &bestsubtype);
      
      SCIPdebugMessage("  %d: %g <%s, idx=%d, lp=%g, [%g(%d),%g(%d)]>:\n", c, rowcoef, SCIPvarGetName(var), probidx, 
         varsolvals[probidx], bestslb, bestslbtype, bestsub, bestsubtype);

      /* mixed integer set cannot be relaxed to 0-1 single node flow set because both simple bounds are -infinity 
       * and infinity, respectively 
       */
      if( SCIPisInfinity(scip, -bestslb) && SCIPisInfinity(scip, bestsub) )
      {
         assert(!(*success));
         goto TERMINATE;
      }
      
      /* get closest lower bound that can be used to define the real variable y'_j in the 0-1 single node flow 
       * relaxation 
       */
      if( !SCIPisInfinity(scip, bestsub) )
      {
         bestlb = bestslb;
         bestlbtype = bestslbtype;

         if( SCIPvarGetType(var) == SCIP_VARTYPE_CONTINUOUS )
         {
            SCIP_Real bestvlb;
            int bestvlbidx;
            
            SCIP_CALL( getClosestVlb(scip, var, bestsub, rowcoef, rowcoefsbinary, varsolvals, assoctransvars, &bestvlb, &bestvlbidx) );
            if( SCIPisGT(scip, bestvlb, bestlb) )
            {
               bestlb = bestvlb;
               bestlbtype = bestvlbidx;
            }
         }
      }
      /* get closest upper bound that can be used to define the real variable y'_j in the 0-1 single node flow 
       * relaxation 
       */
      if( !SCIPisInfinity(scip, -bestslb) )
      {
         bestub = bestsub;
         bestubtype = bestsubtype;
         
         if( SCIPvarGetType(var) == SCIP_VARTYPE_CONTINUOUS )
         {
            SCIP_Real bestvub;
            int bestvubidx;
         
            SCIP_CALL( getClosestVub(scip, var, bestslb, rowcoef, rowcoefsbinary, varsolvals, assoctransvars, &bestvub, &bestvubidx) );
            if( SCIPisLT(scip, bestvub, bestub) )
            {
               bestub = bestvub;
               bestubtype = bestvubidx;
            }
         }
      }
      SCIPdebugMessage("        bestlb=%g(%d), bestub=%g(%d)\n", bestlb, bestlbtype, bestub, bestubtype);
      
      /* mixed integer set cannot be relaxed to 0-1 single node flow set because there are no suitable bounds 
       * to define the transformed variable y'_j 
       */       
      if( SCIPisInfinity(scip, -bestlb) && SCIPisInfinity(scip, bestub) )
      {
         assert(!(*success));
         goto TERMINATE;
      }
      
      /* select best upper bound if it is closer to the LP value of y_j and best lower bound otherwise and use this bound 
       * to define the real variable y'_j with 0 <= y'_j <= u'_j x_j in the 0-1 single node flow relaxation; 
       * prefer variable bounds 
       */
      if( SCIPisEQ(scip, varsolvals[probidx], (1.0 - BOUNDSWITCH) * bestlb + BOUNDSWITCH * bestub) && bestlbtype >= 0 )
         uselb = TRUE;
      else if( SCIPisEQ(scip, varsolvals[probidx], (1.0 - BOUNDSWITCH) * bestlb + BOUNDSWITCH * bestub) 
         && bestubtype >= 0 )
         uselb = FALSE;
      else if( SCIPisLE(scip, varsolvals[probidx], (1.0 - BOUNDSWITCH) * bestlb + BOUNDSWITCH * bestub) )
         uselb = TRUE;
      else
      {
         assert(SCIPisGT(scip, varsolvals[probidx], (1.0 - BOUNDSWITCH) * bestlb + BOUNDSWITCH * bestub));
         uselb = FALSE;
      }
      if( uselb )
      {
         /* use bestlb to define y'_j */

         assert(!SCIPisInfinity(scip, bestsub));
         assert(!SCIPisInfinity(scip, - bestlb));
         assert(bestsubtype == -1 || bestsubtype == -2);
         assert(bestlbtype > -3 && bestlbtype < SCIPvarGetNVlbs(var));

         /* store for y_j that bestlb is the bound used to define y'_j and that y'_j is the associated real variable 
          * in the relaxed set 
          */ 
         boundsfortrans[probidx] = bestlbtype;
         boundtypesfortrans[probidx] = SCIP_BOUNDTYPE_LOWER;
         assoctransvars[probidx] = *ntransvars;
   
         if( bestlbtype < 0 )
         {
            SCIP_Real val;
            SCIP_Real contsolval;

            /* use simple lower bound in bestlb = l_j <= y_j <= u_j = bestsub to define 
             *   y'_j = - a_j ( y_j - u_j ) with 0 <= y'_j <=   a_j ( u_j - l_j ) x_j and x_j = 1    if a_j > 0
             *   y'_j =   a_j ( y_j - u_j ) with 0 <= y'_j <= - a_j ( u_j - l_j ) x_j and x_j = 1    if a_j < 0,
             * put j into the set 
             *   N2   if a_j > 0
             *   N1   if a_j < 0
             * and update the right hand side of the constraint in the relaxation
             *   rhs = rhs - a_j u_j
             */
            val = rowcoef * ( bestsub - bestlb );
            contsolval = rowcoef * ( varsolvals[probidx] - bestsub );
            if( SCIPisPositive(scip, rowcoef) )
            {
               transvarcoefs[*ntransvars] = - 1;
               transvarvubcoefs[*ntransvars] = val;
               transbinvarsolvals[*ntransvars] = 1.0;
               transcontvarsolvals[*ntransvars] = - contsolval;
            }
            else
            {
               assert(SCIPisNegative(scip, rowcoef));
               transvarcoefs[*ntransvars] = 1;
               transvarvubcoefs[*ntransvars] = - val;
               transbinvarsolvals[*ntransvars] = 1.0;
               transcontvarsolvals[*ntransvars] = contsolval;
            }
            (*transrhs) -= (rowcoef * bestsub);

            SCIPdebugMessage("    --> bestlb used for trans: ... %s y'_%d + ..., y'_%d <= %g x_%d (=1), rhs=%g-(%g*%g)=%g\n", 
               transvarcoefs[*ntransvars] == 1 ? "+" : "-", *ntransvars, *ntransvars, transvarvubcoefs[*ntransvars], 
               *ntransvars, *transrhs + (rowcoef * bestsub), rowcoef, bestsub, *transrhs);
         }
         else
         {
            SCIP_Real rowcoefbinary;
            SCIP_Real varsolvalbinary;
            SCIP_Real val;
            SCIP_Real contsolval;
            SCIP_VAR** vlbvars = SCIPvarGetVlbVars(var);
            SCIP_Real* vlbconsts = SCIPvarGetVlbConstants(var);
            SCIP_Real* vlbcoefs = SCIPvarGetVlbCoefs(var);
            
            /* use variable lower bound in bestlb = l~_j x_j + d_j <= y_j <= u_j = bestsub to define 
             *   y'_j = - ( a_j ( y_j - d_j ) + c_j x_j ) with 0 <= y'_j <= - ( a_j l~_j + c_j ) x_j    if a_j > 0
             *   y'_j =     a_j ( y_j - d_j ) + c_j x_j   with 0 <= y'_j <=   ( a_j l~_j + c_j ) x_j    if a_j < 0,
             * where c_j is the coefficient of x_j in the row, put j into the set 
             *   N2   if a_j > 0
             *   N1   if a_j < 0
             * and update the right hand side of the constraint in the relaxation
             *   rhs = rhs - a_j d_j
             */
            assert(SCIPvarIsBinary(vlbvars[bestlbtype]));

            rowcoefbinary = rowcoefsbinary[SCIPvarGetProbindex(vlbvars[bestlbtype])];
            varsolvalbinary = varsolvals[SCIPvarGetProbindex(vlbvars[bestlbtype])];

            val = (rowcoef * vlbcoefs[bestlbtype]) + rowcoefbinary;
            contsolval = (rowcoef * (varsolvals[probidx] - vlbconsts[bestlbtype])) + (rowcoefbinary * varsolvalbinary);

            if( SCIPisPositive(scip, rowcoef) )
            {
               transvarcoefs[*ntransvars] = - 1;
               transvarvubcoefs[*ntransvars] = - val;
               transbinvarsolvals[*ntransvars] = varsolvalbinary;
               transcontvarsolvals[*ntransvars] = - contsolval;
            }
            else
            {
               assert(SCIPisNegative(scip, rowcoef));
               transvarcoefs[*ntransvars] = 1;
               transvarvubcoefs[*ntransvars] = val;
               transbinvarsolvals[*ntransvars] = varsolvalbinary;
               transcontvarsolvals[*ntransvars] = contsolval;
            }
            (*transrhs) -= (rowcoef * vlbconsts[bestlbtype]); 

            /* store for x_j that y'_j is the associated real variable in the 0-1 single node flow relaxation */
            assoctransvars[SCIPvarGetProbindex(vlbvars[bestlbtype])] = *ntransvars;

            SCIPdebugMessage("    --> bestlb used for trans: ... %s y'_%d + ..., y'_%d <= %g x_%d (=%s), rhs=%g-(%g*%g)=%g\n", 
               transvarcoefs[*ntransvars] == 1 ? "+" : "-", *ntransvars, *ntransvars, transvarvubcoefs[*ntransvars], 
               *ntransvars, SCIPvarGetName(vlbvars[bestlbtype]), *transrhs + (rowcoef * vlbconsts[bestlbtype]), rowcoef, 
               vlbconsts[bestlbtype], *transrhs );
         }
      }
      else
      {
         /* use bestub to define y'_j */

         assert(!SCIPisInfinity(scip, bestub));
         assert(!SCIPisInfinity(scip, - bestslb));
         assert(bestslbtype == -1 || bestslbtype == -2);
         assert(bestubtype > -3 && bestubtype < SCIPvarGetNVubs(var));

         /* store for y_j that bestub is the bound used to define y'_j and that y'_j is the associated real variable 
          * in the relaxed set 
          */ 
         boundsfortrans[probidx] = bestubtype;
         boundtypesfortrans[probidx] = SCIP_BOUNDTYPE_UPPER;
         assoctransvars[probidx] = *ntransvars;
   
         if( bestubtype < 0 )
         {
            SCIP_Real val;
            SCIP_Real contsolval;

            /* use simple upper bound in bestslb = l_j <= y_j <= u_j = bestub to define 
             *   y'_j =   a_j ( y_j - l_j ) with 0 <= y'_j <=   a_j ( u_j - l_j ) x_j and x_j = 1    if a_j > 0
             *   y'_j = - a_j ( y_j - l_j ) with 0 <= y'_j <= - a_j ( u_j - l_j ) x_j and x_j = 1    if a_j < 0,
             * put j into the set 
             *   N1   if a_j > 0
             *   N2   if a_j < 0
             * and update the right hand side of the constraint in the relaxation
             *   rhs = rhs - a_j l_j
             */
            val = rowcoef * ( bestub - bestslb );
            contsolval = rowcoef * ( varsolvals[probidx] - bestslb );
            if( SCIPisPositive(scip, rowcoef) )
            {
               transvarcoefs[*ntransvars] = 1;
               transvarvubcoefs[*ntransvars] = val;
               transbinvarsolvals[*ntransvars] = 1.0;
               transcontvarsolvals[*ntransvars] = contsolval;
            }
            else
            {
               assert(SCIPisNegative(scip, rowcoef));
               transvarcoefs[*ntransvars] = - 1;
               transvarvubcoefs[*ntransvars] = - val;
               transbinvarsolvals[*ntransvars] = 1.0;
               transcontvarsolvals[*ntransvars] = - contsolval;
            }
            (*transrhs) -= (rowcoef * bestslb);

            SCIPdebugMessage("    --> bestub used for trans: ... %s y'_%d + ..., Y'_%d <= %g x_%d (=1), rhs=%g-(%g*%g)=%g\n", 
               transvarcoefs[*ntransvars] == 1 ? "+" : "-", *ntransvars, *ntransvars, transvarvubcoefs[*ntransvars], 
               *ntransvars, *transrhs + (rowcoef * bestslb), rowcoef, bestslb, *transrhs);
         }
         else
         {
            SCIP_Real rowcoefbinary;
            SCIP_Real varsolvalbinary;
            SCIP_Real val;
            SCIP_Real contsolval;

            SCIP_VAR** vubvars = SCIPvarGetVubVars(var);
            SCIP_Real* vubconsts = SCIPvarGetVubConstants(var);
            SCIP_Real* vubcoefs = SCIPvarGetVubCoefs(var);

            /* use variable upper bound in bestslb = l_j <= y_j <= u~_j x_j + d_j = bestub to define 
             *   y'_j =     a_j ( y_j - d_j ) + c_j x_j   with 0 <= y'_j <=   ( a_j u~_j + c_j ) x_j    if a_j > 0
             *   y'_j = - ( a_j ( y_j - d_j ) + c_j x_j ) with 0 <= y'_j <= - ( a_j u~_j + c_j ) x_j    if a_j < 0,
             * where c_j is the coefficient of x_j in the row, put j into the set 
             *   N1   if a_j > 0
             *   N2   if a_j < 0
             * and update the right hand side of the constraint in the relaxation
             *   rhs = rhs - a_j d_j
             */
            assert(SCIPvarIsBinary(vubvars[bestubtype]));

            rowcoefbinary = rowcoefsbinary[SCIPvarGetProbindex(vubvars[bestubtype])];
            varsolvalbinary = varsolvals[SCIPvarGetProbindex(vubvars[bestubtype])];

            val = ( rowcoef * vubcoefs[bestubtype] ) + rowcoefbinary;
            contsolval = (rowcoef * (varsolvals[probidx] - vubconsts[bestubtype])) + (rowcoefbinary * varsolvalbinary);

            if( SCIPisPositive(scip, rowcoef) )
            {
               transvarcoefs[*ntransvars] = 1;
               transvarvubcoefs[*ntransvars] = val;
               transbinvarsolvals[*ntransvars] = varsolvalbinary;
               transcontvarsolvals[*ntransvars] = contsolval;
            }
            else
            {
               assert(SCIPisNegative(scip, rowcoef));
               transvarcoefs[*ntransvars] = - 1;
               transvarvubcoefs[*ntransvars] = - val;
               transbinvarsolvals[*ntransvars] = varsolvalbinary;
               transcontvarsolvals[*ntransvars] = - contsolval;
            }
            (*transrhs) -= (rowcoef * vubconsts[bestubtype]); 

            /* store for x_j that y'_j is the associated real variable in the 0-1 single node flow relaxation */
            assoctransvars[SCIPvarGetProbindex(vubvars[bestubtype])] = *ntransvars;

            SCIPdebugMessage("    --> bestub used for trans: ... %s y'_%d + ..., y'_%d <= %g x_%d (=%s), rhs=%g-(%g*%g)=%g\n", 
               transvarcoefs[*ntransvars] == 1 ? "+" : "-", *ntransvars, *ntransvars, transvarvubcoefs[*ntransvars], 
               *ntransvars, SCIPvarGetName(vubvars[bestubtype]), *transrhs + (rowcoef * vubconsts[bestubtype]), rowcoef, 
               vubconsts[bestubtype], *transrhs);
         }
      }

      /* relaxing the mixed integer set to a 0-1 single node flow set was not successful because coefficient of y_j and
       * the bounds selected for the transformation together result in an infinite variable upper bound in the 0-1 single
       * node flow set; this can be caused by huge but finite values for the bounds or the coefficient
       */
      if( SCIPisInfinity(scip, transvarvubcoefs[*ntransvars]) )
      {
         assert(!(*success));
         goto TERMINATE;
      }

      assert(boundsfortrans[probidx] > -3);
      assert(assoctransvars[probidx] >= 0 && assoctransvars[probidx] == (*ntransvars));
      assert(transvarcoefs[*ntransvars] == 1 || transvarcoefs[*ntransvars] == - 1 );
      assert(SCIPisFeasGE(scip, transbinvarsolvals[*ntransvars], 0.0) 
         && SCIPisFeasLE(scip, transbinvarsolvals[*ntransvars], 1.0));
      assert(SCIPisFeasGE(scip, transvarvubcoefs[*ntransvars], 0.0) 
         && !SCIPisInfinity(scip, transvarvubcoefs[*ntransvars]));

      /* updates number of variables in transformed problem */
      (*ntransvars)++;
   }
   assert(*ntransvars == nnonzcolsnonbinary);

   SCIPdebugMessage("transformation for BINARY variables (nbinvars=%d):\n", nnonzcolsbinary);

   /* binary variables not involved in used variable bounds imposed on non-binary variable */
   for( c = 0; c < nnonzcolsbinary; c++ )
   {
      SCIP_VAR* var;
      SCIP_Real rowcoef; 
      int probidx;
      SCIP_Real val;
      SCIP_Real contsolval;
      
      var = SCIPcolGetVar(nonzcols[nonzcolsbinary[c]]);
      probidx = SCIPvarGetProbindex(var);
      rowcoef = rowweight * scale * nonzcoefs[nonzcolsbinary[c]];

      assert(rowcoefsbinary[probidx] == rowcoef); /*lint !e777*/
      assert(!SCIPisZero(scip, rowcoef));

      SCIPdebugMessage("  %d: %g <%s, idx=%d, lp=%g, [%g, %g]>:\n", c, rowcoef, SCIPvarGetName(var), probidx, varsolvals[probidx], 
         SCIPvarGetLbGlobal(var), SCIPvarGetUbGlobal(var));
      
      /* x_j has already been handled in connection with a non-binary variable */ 
      if( assoctransvars[probidx] > -1 )
      {
         assert(assoctransvars[probidx] >= 0 && assoctransvars[probidx] <= nnonzcolsnonbinary);
         assert(boundsfortrans[probidx] == -3);
         SCIPdebugMessage("   --> already handled\n");
         continue;
      }
      
      assert(assoctransvars[probidx] == -1);
      assert(boundsfortrans[probidx] == -3);

      /* store for x_j that y'_j is the associated real variable in the 0-1 single node flow relaxation */
      assoctransvars[probidx] = *ntransvars;
      
      /* define
       *    y'_j =   c_j x_j with 0 <= y'_j <=   c_j x_j    if c_j > 0   
       *    y'_j = - c_j x_j with 0 <= y'_j <= - c_j x_j    if c_j < 0,   
       * where c_j is the coefficient of x_j in the row and put j into the set
       *    N1   if c_j > 0
       *    N2   if c_j < 0.
       */
      val = rowcoef;
      contsolval = rowcoef * varsolvals[probidx];
      if( SCIPisPositive(scip, rowcoef) )
      {
         transvarcoefs[*ntransvars] = 1;
         transvarvubcoefs[*ntransvars] = val;
         transbinvarsolvals[*ntransvars] = varsolvals[probidx];
         transcontvarsolvals[*ntransvars] = contsolval;
      }
      else
      {
         assert(SCIPisNegative(scip, rowcoef));
         transvarcoefs[*ntransvars] = - 1;
         transvarvubcoefs[*ntransvars] = - val;
         transbinvarsolvals[*ntransvars] = varsolvals[probidx];
         transcontvarsolvals[*ntransvars] = - contsolval;
      }
      assert(assoctransvars[probidx] >= 0 && assoctransvars[probidx] == (*ntransvars));
      assert(transvarcoefs[*ntransvars] == 1 || transvarcoefs[*ntransvars] == - 1 );
      assert(SCIPisFeasGE(scip, transbinvarsolvals[*ntransvars], 0.0) 
         && SCIPisFeasLE(scip, transbinvarsolvals[*ntransvars], 1.0));
      assert(SCIPisFeasGE(scip, transvarvubcoefs[*ntransvars], 0.0) 
         && !SCIPisInfinity(scip, transvarvubcoefs[*ntransvars]));
         
      SCIPdebugMessage("   --> ... %s y'_%d + ..., y'_%d <= %g x_%d (=%s))\n", 
         transvarcoefs[*ntransvars] == 1 ? "+" : "-", *ntransvars, *ntransvars, 
         transvarvubcoefs[*ntransvars], *ntransvars, SCIPvarGetName(var) );

      /* updates number of variables in transformed problem */
      (*ntransvars)++;
   }
   assert(*ntransvars >= nnonzcolsnonbinary && *ntransvars <= nnonzcols);

   /* construction was successful */
   *success = TRUE;

#ifdef SCIP_DEBUG
   SCIPdebugMessage("constraint in constructed 0-1 single node flow relaxation: ");
   for( c = 0; c < *ntransvars; c++ )
   {   
      SCIPdebugPrintf("%s y'_%d ", transvarcoefs[c] == 1 ? "+" : "-", c);
   }
   SCIPdebugPrintf("<= %g\n", *transrhs);
#endif

 TERMINATE:
   /* free data structures */
   SCIPfreeBufferArray(scip, &rowcoefsbinary);
   SCIPfreeBufferArray(scip, &nonzcolsnonbinary);
   SCIPfreeBufferArray(scip, &nonzcolsbinary);

   return SCIP_OKAY;
}

/** solve knapsack problem in maximization form with "<" constraint approximately by greedy; if needed, one can provide 
 *  arrays to store all selected items and all not selected items
 */
static
SCIP_RETCODE SCIPsolveKnapsackApproximatelyLT(
   SCIP*                 scip,               /**< SCIP data structure */
   int                   nitems,             /**< number of available items */
   SCIP_Real*            weights,            /**< item weights */
   SCIP_Real*            profits,            /**< item profits */
   SCIP_Real             capacity,           /**< capacity of knapsack */
   int*                  items,              /**< item numbers */
   int*                  solitems,           /**< array to store items in solution, or NULL */
   int*                  nonsolitems,        /**< array to store items not in solution, or NULL */
   int*                  nsolitems,          /**< pointer to store number of items in solution, or NULL */
   int*                  nnonsolitems,       /**< pointer to store number of items not in solution, or NULL */
   SCIP_Real*            solval              /**< pointer to store optimal solution value, or NULL */
   ) 
{
   SCIP_Real* tempsort;
   SCIP_Real solitemsweight;
   int j;
   int i;
   
   assert(weights != NULL);
   assert(profits != NULL);
   assert(SCIPisFeasGE(scip, capacity, 0.0));
   assert(!SCIPisInfinity(scip, capacity));
   assert(items != NULL);
   assert(nitems >= 0);

   if( solitems != NULL )
   {
      *nsolitems = 0;
      *nnonsolitems = 0;
   }
   if( solval != NULL )
      *solval = 0.0;

   /* allocate memory for temporary array used for sorting; array should contain profits divided by corresponding weights (p_1 / w_1 ... p_n / w_n )*/
   SCIP_CALL( SCIPallocBufferArray(scip, &tempsort, nitems) );
   /* initialize temporary array */ 
   for( i = nitems - 1; i >= 0; --i )
   {
      tempsort[i] = profits[i] / weights [i];
   }

   /* sort tempsort, items, weights and profits such that p_1 / w_1 >= p_2 / w_2 >= ... >= p_n / w_n */
   SCIPsortDownRealRealRealInt ( tempsort, weights, profits, items, nitems);

   /* free temporary array */
   SCIPfreeBufferArray(scip, &tempsort);
   
   /* select items as long as they fit into the knapsack */
   solitemsweight = 0.0;
   for( j = 0; j < nitems && SCIPisFeasLT(scip, solitemsweight + weights[j], capacity); j++ )
   {
      if( solitems != NULL )
      {
         solitems[*nsolitems] = items[j];
         (*nsolitems)++;
      }
      if( solval != NULL )
         (*solval) += profits[j];
      solitemsweight += weights[j];
   }
   for( ; j < nitems && solitems != NULL; j++ )
   {
      nonsolitems[*nnonsolitems] = items[j];
      (*nnonsolitems)++;
   }
   
   return SCIP_OKAY;
}

/** checks, whether the given scalar scales the given value to an integral number with error in the given bounds */
static
SCIP_Bool isIntegralScalar(
   SCIP_Real             val,                /**< value that should be scaled to an integral value */
   SCIP_Real             scalar,             /**< scalar that should be tried */
   SCIP_Real             mindelta,           /**< minimal relative allowed difference of scaled coefficient s*c and integral i */
   SCIP_Real             maxdelta            /**< maximal relative allowed difference of scaled coefficient s*c and integral i */
   )
{
   SCIP_Real sval;
   SCIP_Real downval;
   SCIP_Real upval;

   assert(mindelta <= 0.0);
   assert(maxdelta >= 0.0);

   sval = val * scalar;
   downval = floor(sval);
   upval = ceil(sval);

   return (SCIPrelDiff(sval, downval) <= maxdelta || SCIPrelDiff(sval, upval) >= mindelta);
}

/** get integral number with error in the bounds which corresponds to given value scaled by a given scalar;
 *  should be used in connection with isIntegralScalar()  
 */
static 
SCIP_Longint getIntegralVal(
   SCIP_Real             val,                /**< value that should be scaled to an integral value */
   SCIP_Real             scalar,             /**< scalar that should be tried */
   SCIP_Real             mindelta,           /**< minimal relative allowed difference of scaled coefficient s*c and integral i */
   SCIP_Real             maxdelta            /**< maximal relative allowed difference of scaled coefficient s*c and integral i */
   )
{
   SCIP_Real sval;
   SCIP_Real upval;
   SCIP_Longint intval;

   assert(mindelta <= 0.0);
   assert(maxdelta >= 0.0);

   sval = val * scalar;
   upval = ceil(sval);

   if( SCIPrelDiff(sval, upval) >= mindelta )
      intval = (SCIP_Longint) upval;
   else
      intval = (SCIP_Longint) (floor(sval));
   
   return intval;
}

/** build the flow cover which corresponds to the given exact or approximate solution of KP^SNF; given unfinished 
 *  flow cover contains variables which have been fixed in advance 
 */  
static
void buildFlowCover(
   SCIP*                 scip,               /**< SCIP data structure */ 
   int*                  coefs,              /**< coefficient of all real variables in N1&N2 */ 
   SCIP_Real*            vubcoefs,           /**< coefficient in vub of all real variables in N1&N2 */
   SCIP_Real             rhs,                /**< right hand side of 0-1 single node flow constraint */ 
   int*                  solitems,           /**< items in knapsack */
   int*                  nonsolitems,        /**< items not in knapsack */
   int                   nsolitems,          /**< number of items in knapsack */
   int                   nnonsolitems,       /**< number of items not in knapsack */
   int*                  nflowcovervars,     /**< pointer to store number of variables in flow cover */
   int*                  nnonflowcovervars,  /**< pointer to store number of variables not in flow cover */
   int*                  flowcoverstatus,    /**< pointer to store whether variable is in flow cover (+1) or not (-1) */ 
   SCIP_Real*            flowcoverweight,    /**< pointer to store weight of flow cover */
   SCIP_Real*            lambda              /**< pointer to store lambda */
)
{
   int j; 

   assert(scip != NULL);
   assert(coefs != NULL);
   assert(vubcoefs != NULL);
   assert(solitems != NULL);
   assert(nonsolitems != NULL);
   assert(nsolitems >= 0);
   assert(nnonsolitems >= 0);
   assert(nflowcovervars != NULL && *nflowcovervars >= 0);
   assert(nnonflowcovervars != NULL && *nnonflowcovervars >= 0);
   assert(flowcoverstatus != NULL);
   assert(flowcoverweight != NULL);
   assert(lambda != NULL);

   /* get flowcover status for each item */
   for( j = 0; j < nsolitems; j++ )
   {
      /* j in N1 with z°_j = 1 => j in N1\C1 */
      if( coefs[solitems[j]] == 1 )
      {
         flowcoverstatus[solitems[j]] = -1;
         (*nnonflowcovervars)++;
      }
      /* j in N2 with z_j = 1 => j in C2 */
      else
      {
         assert(coefs[solitems[j]] == -1);
         flowcoverstatus[solitems[j]] = 1;
         (*nflowcovervars)++;
         (*flowcoverweight) -= vubcoefs[solitems[j]];
      }
   }
   for( j = 0; j < nnonsolitems; j++ )
   {
      /* j in N1 with z°_j = 0 => j in C1 */
      if( coefs[nonsolitems[j]] == 1 )
      {
         flowcoverstatus[nonsolitems[j]] = 1;
         (*nflowcovervars)++;
         (*flowcoverweight) += vubcoefs[nonsolitems[j]];
      }
      /* j in N2 with z_j = 0 => j in N2\C2 */
      else
      {
         assert(coefs[nonsolitems[j]] == -1);
         flowcoverstatus[nonsolitems[j]] = -1;
         (*nnonflowcovervars)++;
      }
   }

   /* get lambda = sum_{j in C1} u_j - sum_{j in C2} u_j - rhs */
   *lambda = (*flowcoverweight) - rhs;
}

/** get a flow cover (C1, C2) for a given 0-1 single node flow set 
 *    {(x,y) in {0,1}^n x R^n : sum_{j in N1} y_j - sum_{j in N2} y_j <= b, 0 <= y_j <= u_j x_j}, 
 *  i.e., get sets C1 subset N1 and C2 subset N2 with sum_{j in C1} u_j - sum_{j in C2} u_j = b + lambda and lambda > 0  
 */
static
SCIP_RETCODE getFlowCover(
   SCIP*                 scip,               /**< SCIP data structure */ 
   int*                  coefs,              /**< coefficient of all real variables in N1&N2 */ 
   SCIP_Real*            solvals,            /**< LP solution value of binary variable in vub of all real vars in N1&N2 */
   SCIP_Real*            vubcoefs,           /**< coefficient in vub of all real variables in N1&N2 */
   int                   nvars,              /**< number of real variables in N1&N2 */
   SCIP_Real             rhs,                /**< right hand side of 0-1 single node flow constraint */ 
   int*                  nflowcovervars,     /**< pointer to store number of variables in flow cover */
   int*                  nnonflowcovervars,  /**< pointer to store number of variables not in flow cover */
   int*                  flowcoverstatus,    /**< pointer to store whether variable is in flow cover (+1) or not (-1) */ 
   SCIP_Real*            lambda,             /**< pointer to store lambda */ 
   SCIP_Bool*            found               /**< pointer to store whether a cover was found */
   )
{
   SCIP_Real* transprofitsint;
   SCIP_Real* transprofitsreal;
   SCIP_Real* transweightsreal;
   SCIP_Longint* transweightsint;
   int* items;
   int* itemsint;
   int* nonsolitems;
   int* solitems;
   SCIP_Real flowcoverweight;
   SCIP_Real flowcoverweightafterfix;
   SCIP_Real n1itemsweight;
   SCIP_Real n2itemsminweight;
   SCIP_Real scalar;
   SCIP_Real transcapacityreal;
#if !defined(NDEBUG) || defined(SCIP_DEBUG)
   SCIP_Bool kpexact;
#endif
   SCIP_Bool scalesuccess;
   SCIP_Bool transweightsrealintegral;
   SCIP_Longint transcapacityint;
   int nflowcovervarsafterfix;
   int nitems;
   int nn1items;
   int nnonflowcovervarsafterfix;
   int nnonsolitems;
   int nsolitems;
   int j;
   
   assert(scip != NULL);
   assert(coefs != NULL);
   assert(solvals != NULL);
   assert(vubcoefs != NULL);
   assert(nvars > 0);
   assert(nflowcovervars != NULL);
   assert(nnonflowcovervars != NULL);
   assert(flowcoverstatus != NULL);
   assert(lambda != NULL);
   assert(found != NULL);

   SCIPdebugMessage("--------------------- get flow cover ----------------------------------------------------\n");

   /* get data structures */
   SCIP_CALL( SCIPallocBufferArray(scip, &items, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &itemsint, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &transprofitsreal, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &transprofitsint, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &transweightsreal, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &transweightsint, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &solitems, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &nonsolitems, nvars) );

   BMSclearMemoryArray(flowcoverstatus, nvars);
   *found = FALSE;
   *nflowcovervars = 0;
   *nnonflowcovervars = 0;
   flowcoverweight = 0.0;
   nflowcovervarsafterfix = 0;
   nnonflowcovervarsafterfix = 0;
   flowcoverweightafterfix = 0.0;
#if !defined(NDEBUG) || defined(SCIP_DEBUG)
   kpexact = FALSE;
#endif

   /* fix some variables in advance according to the following fixing strategy
    *   put j into N1\C1,          if j in N1 and x*_j = 0, 
    *   put j into C1,             if j in N1 and x*_j = 1, 
    *   put j into C2,             if j in N2 and x*_j = 1, 
    *   put j into N2\C2,          if j in N2 and x*_j = 0 
    * and get the set of the remaining variables
    */
   SCIPdebugMessage("0. Fix some variables in advance:\n");
   nitems = 0;
   nn1items = 0;
   n1itemsweight = 0.0;
   n2itemsminweight = SCIP_REAL_MAX;
   for( j = 0; j < nvars; j++ )
   {
      assert(coefs[j] == 1 || coefs[j] == -1);
      assert(SCIPisFeasGE(scip, solvals[j], 0.0) && SCIPisFeasLE(scip, solvals[j], 1.0));
      assert(SCIPisFeasGE(scip, vubcoefs[j], 0.0));
         
      /* if u_j = 0, put j into N1\C1 and N2\C2, respectively */
      if( SCIPisFeasZero(scip, vubcoefs[j]) )
      {
         flowcoverstatus[j] = -1;
         (*nnonflowcovervars)++;
         continue;
      }
         
      /* x*_j is fractional */
      if( !SCIPisFeasIntegral(scip, solvals[j]) ) 
      {
         items[nitems] = j;
         nitems++;
         if( coefs[j] == 1 )
         {
            n1itemsweight += vubcoefs[j];
            nn1items++;
         }
         else
            n2itemsminweight = MIN(n2itemsminweight, vubcoefs[j]);
      }
      /* j is in N1 and x*_j = 0 */
      else if( coefs[j] == 1 && solvals[j] < 0.5 )
      {
         flowcoverstatus[j] = -1;
         (*nnonflowcovervars)++;
         SCIPdebugMessage("     <%d>: in N1-C1\n", j);
      }
      /* j is in N1 and x*_j = 1 */
      else if( coefs[j] == 1 && solvals[j] > 0.5 )
      {
         flowcoverstatus[j] = 1;
         (*nflowcovervars)++;
         flowcoverweight += vubcoefs[j];
         SCIPdebugMessage("     <%d>: in C1\n", j);
      }
      /* j is in N2 and x*_j = 1 */
      else if( coefs[j] == -1 && solvals[j] > 0.5 )
      {
         flowcoverstatus[j] = 1;
         (*nflowcovervars)++;
         flowcoverweight -= vubcoefs[j];
         SCIPdebugMessage("     <%d>: in C2\n", j);
      } 
      /* j is in N2 and x*_j = 0 */
      else 
      {
         assert(coefs[j] == -1 && solvals[j] < 0.5);
         flowcoverstatus[j] = -1;
         (*nnonflowcovervars)++;
         SCIPdebugMessage("     <%d>: in N2-C2\n", j);
      }
   }
   assert((*nflowcovervars) + (*nnonflowcovervars) + nitems == nvars);
   assert(nn1items >= 0);
   
   /* to find a flow cover, transform the following knapsack problem
    *
    * (KP^SNF)      max sum_{j in N1} ( x*_j - 1 ) z_j + sum_{j in N2} x*_j z_j
    *                   sum_{j in N1}          u_j z_j - sum_{j in N2} u_j  z_j > b
    *                                         z_j in {0,1} for all j in N1 & N2
    *
    * 1. to a knapsack problem in maximization form, such that all variables in the knapsack constraint have 
    *    positive weights and the constraint is a "<" constraint, by complementing all variables in N1
    *     
    *    (KP^SNF_rat)  max sum_{j in N1} ( 1 - x*_j ) z°_j + sum_{j in N2} x*_j z_j
    *                      sum_{j in N1}          u_j z°_j + sum_{j in N2} u_j  z_j < - b + sum_{j in N1} u_j  
    *                                                 z°_j in {0,1} for all j in N1 
    *                                                  z_j in {0,1} for all j in N2, 
    *    and solve it approximately under consideration of the fixing, 
    * or 
    * 2. to a knapsack problem in maximization form, such that all variables in the knapsack constraint have 
    *    positive integer weights and the constraint is a "<=" constraint, by complementing all variables in N1
    *    and multiplying the constraint by a suitable scalar C
    *
    *    (KP^SNF_int)  max sum_{j in N1} ( 1 - x*_j ) z°_j + sum_{j in N2} x*_j z_j
    *                      sum_{j in N1}        C u_j z°_j + sum_{j in N2} C u_j  z_j <= c 
    *                                                   z°_j in {0,1} for all j in N1 
    *                                                    z_j in {0,1} for all j in N2, 
    *    where 
    *      c = floor[ C (- b + sum_{j in N1} u_j ) ]      if frac[ C (- b + sum_{j in N1} u_j ) ] > 0
    *      c =        C (- b + sum_{j in N1} u_j )   - 1  if frac[ C (- b + sum_{j in N1} u_j ) ] = 0
    *    and solve it exactly under consideration of the fixing.
    */
   SCIPdebugMessage("1. Transform KP^SNF to KP^SNF_rat:\n");
   /* get weight and profit of variables in KP^SNF_rat and check, whether all weights are already integral */
   transweightsrealintegral = TRUE;
   for( j = 0; j < nitems; j++ )
   {
      transweightsreal[j] = vubcoefs[items[j]];

      if( !isIntegralScalar(transweightsreal[j], 1.0, -MINDELTA, MAXDELTA) )
         transweightsrealintegral = FALSE;

      if( coefs[items[j]] == 1 )
      {         
         transprofitsreal[j] = 1.0 - solvals[items[j]];
         SCIPdebugMessage("     <%d>: j in N1:   w_%d = %g, p_%d = %g %s\n", items[j], items[j], transweightsreal[j], 
            items[j], transprofitsreal[j], SCIPisIntegral(scip, transweightsreal[j]) ? "" : "  ----> NOT integral");
      }
      else
      {
         transprofitsreal[j] = solvals[items[j]];
         SCIPdebugMessage("     <%d>: j in N2:   w_%d = %g, p_%d = %g %s\n", items[j], items[j], transweightsreal[j], 
            items[j], transprofitsreal[j], SCIPisIntegral(scip, transweightsreal[j]) ? "" : "  ----> NOT integral");
      }
   }
   /* get capacity of knapsack constraint in KP^SNF_rat */
   transcapacityreal = - rhs + flowcoverweight + n1itemsweight;
   SCIPdebugMessage("     transcapacity = -rhs(%g) + flowcoverweight(%g) + n1itemsweight(%g) = %g\n", 
      rhs, flowcoverweight, n1itemsweight, transcapacityreal);

   /* there exists no flow cover if the capacity of knapsack constraint in KP^SNF_rat after fixing 
    * is less than or equal to zero 
    */ 
   if( SCIPisFeasLE(scip, transcapacityreal/10, 0.0) )
   {
      assert(!(*found));
      goto TERMINATE;
   }

   /* KP^SNF_rat has been solved by fixing some variables in advance */
   assert(nitems >= 0);
   if( nitems == 0)
   {
      /* get lambda = sum_{j in C1} u_j - sum_{j in C2} u_j - rhs */
      *lambda = flowcoverweight - rhs;
      *found = TRUE;
      goto TERMINATE;
   }
   
   /* Use the following strategy 
    *   solve KP^SNF_int exactly,          if a suitable factor C is found and (nitems*capacity) <= MAXDYNPROGSPACE, 
    *   solve KP^SNF_rat approximately,    otherwise 
    */ 
   
   /* find a scaling factor C */
   if( transweightsrealintegral )
   {
      /* weights are already integral */
      scalar = 1.0;
      scalesuccess = TRUE;
   }
   else
   {
      scalesuccess = FALSE;
      SCIP_CALL( SCIPcalcIntegralScalar(transweightsreal, nitems, -MINDELTA, MAXDELTA, MAXDNOM, MAXSCALE, &scalar, 
            &scalesuccess) );
   }

   /* initialize number of (non-)solution items, should be changed to a nonnegative number in all possible paths below */
   nsolitems = -1;
   nnonsolitems = -1;

   /* suitable factor C was found*/
   if( scalesuccess )
   {
      SCIP_Real tmp1;
      SCIP_Real tmp2;
         
      /* transform KP^SNF to KP^SNF_int */
      for( j = 0; j < nitems; ++j )
      {
         transweightsint[j] = getIntegralVal(transweightsreal[j], scalar, -MINDELTA, MAXDELTA);
         transprofitsint[j] = transprofitsreal[j];
         itemsint[j] = items[j];
      }
      if( isIntegralScalar(transcapacityreal, scalar, -MINDELTA, MAXDELTA) )
      {
         transcapacityint = getIntegralVal(transcapacityreal, scalar, -MINDELTA, MAXDELTA);
         transcapacityint -= 1;
      }
      else
         transcapacityint = (SCIP_Longint) (transcapacityreal * scalar);
      nflowcovervarsafterfix = *nflowcovervars;
      nnonflowcovervarsafterfix = *nnonflowcovervars;
      flowcoverweightafterfix = flowcoverweight;

      tmp1 = (SCIP_Real) (nitems + 1);
      tmp2 = (SCIP_Real) ((transcapacityint) + 1);
      if( transcapacityint * nitems <= MAXDYNPROGSPACE && tmp1 * tmp2 <= INT_MAX / 8.0)
      {
         SCIP_Bool success;

         /* solve KP^SNF_int by dynamic programming */
         SCIP_CALL(SCIPsolveKnapsackExactly(scip, nitems, transweightsint, transprofitsint, transcapacityint, 
               itemsint, solitems, nonsolitems, &nsolitems, &nnonsolitems, NULL, &success));

         if( !success )
         {
            /* solve KP^SNF_rat approximately */
            SCIP_CALL(SCIPsolveKnapsackApproximatelyLT(scip, nitems, transweightsreal, transprofitsreal, 
                  transcapacityreal, items, solitems, nonsolitems, &nsolitems, &nnonsolitems, NULL));
         }
#if !defined(NDEBUG) || defined(SCIP_DEBUG)
         else
            kpexact = TRUE;
#endif
      }
      else
      {
         /* solve KP^SNF_rat approximately */
         SCIP_CALL(SCIPsolveKnapsackApproximatelyLT(scip, nitems, transweightsreal, transprofitsreal, transcapacityreal,
               items, solitems, nonsolitems, &nsolitems, &nnonsolitems, NULL));
         assert(!kpexact);
      }
   }
   else
   {
      /* solve KP^SNF_rat approximately */
      SCIP_CALL(SCIPsolveKnapsackApproximatelyLT(scip, nitems, transweightsreal, transprofitsreal, transcapacityreal, 
            items, solitems, nonsolitems, &nsolitems, &nnonsolitems, NULL));
      assert(!kpexact);
   }

   assert(nsolitems != -1);
   assert(nnonsolitems != -1);

   /* build the flow cover from the solution of KP^SNF_rat and KP^SNF_int, respectively and the fixing */
   assert(*nflowcovervars + *nnonflowcovervars + nsolitems + nnonsolitems == nvars);
   buildFlowCover(scip, coefs, vubcoefs, rhs, solitems, nonsolitems, nsolitems, nnonsolitems, nflowcovervars,
      nnonflowcovervars, flowcoverstatus, &flowcoverweight, lambda);
   assert(*nflowcovervars + *nnonflowcovervars == nvars);
 
   /* if the found structure is not a flow cover, because of scaling, solve KP^SNF_rat approximately */ 
   if( SCIPisFeasLE(scip, *lambda, 0.0) )
   {
      assert(kpexact);

      /* solve KP^SNF_rat approximately */
      SCIP_CALL(SCIPsolveKnapsackApproximatelyLT(scip, nitems, transweightsreal, transprofitsreal, transcapacityreal, 
            items, solitems, nonsolitems, &nsolitems, &nnonsolitems, NULL));
#ifdef SCIP_DEBUG /* this time only for SCIP_DEBUG, because only then, the variable is used again  */
      kpexact = FALSE;
#endif
      
      /* build the flow cover from the solution of KP^SNF_rat and the fixing */
      *nflowcovervars = nflowcovervarsafterfix;
      *nnonflowcovervars = nnonflowcovervarsafterfix;
      flowcoverweight = flowcoverweightafterfix;
   
      assert(*nflowcovervars + *nnonflowcovervars + nsolitems + nnonsolitems == nvars);
      buildFlowCover(scip, coefs, vubcoefs, rhs, solitems, nonsolitems, nsolitems, nnonsolitems, nflowcovervars,
         nnonflowcovervars, flowcoverstatus, &flowcoverweight, lambda);
      assert(*nflowcovervars + *nnonflowcovervars == nvars);
   }
   *found = TRUE;

 TERMINATE:
   assert((!*found) || SCIPisFeasGT(scip, *lambda, 0.0));  
#ifdef SCIP_DEBUG
   if( *found )
   {
      SCIPdebugMessage("2. %s solution:\n", kpexact ? "exact" : "approximate");
      for( j = 0; j < nvars; j++ )
      {
         if( coefs[j] == 1 && flowcoverstatus[j] == 1 )
         {
            SCIPdebugMessage("     C1: + y_%d [u_%d = %g]\n", j, j, vubcoefs[j]);
         }
         else if( coefs[j] == -1 && flowcoverstatus[j] == 1 )
         {
            SCIPdebugMessage("     C2: - y_%d [u_%d = %g]\n", j, j, vubcoefs[j]);
         }
      }
      SCIPdebugMessage("     flowcoverweight(%g) = rhs(%g) + lambda(%g)\n", flowcoverweight, rhs, *lambda);
   }
#endif
      
   /* free data structures */
   SCIPfreeBufferArray(scip, &nonsolitems);
   SCIPfreeBufferArray(scip, &solitems);
   SCIPfreeBufferArray(scip, &transweightsint);
   SCIPfreeBufferArray(scip, &transweightsreal);
   SCIPfreeBufferArray(scip, &transprofitsint);
   SCIPfreeBufferArray(scip, &transprofitsreal);
   SCIPfreeBufferArray(scip, &itemsint);
   SCIPfreeBufferArray(scip, &items);

   return SCIP_OKAY;
}

/** for a given flow cover and a given value of delta, choose L1 subset N1 \ C1 and L2 subset N2 \ C2 by comparison such that 
 *  the violation of the resulting c-MIRFCI is maximized.
 */
static
void getL1L2(
   SCIP*                 scip,               /**< SCIP data structure */
   int                   ntransvars,         /**< number of continuous variables in N1 & N2 */
   int*                  transvarcoefs,      /**< coefficient of all continuous variables in N1 & N2 */ 
   SCIP_Real*            transbinvarsolvals, /**< LP solution value of bin var in vub of all continuous vars in N1 & N2 */
   SCIP_Real*            transcontvarsolvals,/**< LP solution value of all continuous vars in N1 & N2 */
   SCIP_Real*            transvarvubcoefs,   /**< coefficient of vub of all continuous variables in N1 & N2 */
   int*                  transvarflowcoverstatus,/**< pointer to store whether non-binary var is in L2 (2) or not (-1 or 1) */ 
   SCIP_Real             delta,              /**< delta */
   SCIP_Real             lambda              /**< lambda */
   )
{
   SCIP_Real fbeta;
   SCIP_Real onedivoneminusfbeta;
   int j;

   assert(scip != NULL);
   assert(ntransvars >= 0);
   assert(transvarcoefs != NULL);
   assert(transbinvarsolvals != NULL);
   assert(transcontvarsolvals != NULL);
   assert(transvarvubcoefs != NULL);
   assert(transvarflowcoverstatus != NULL);
   assert(SCIPisGT(scip, delta, lambda));
   assert(SCIPisFeasGT(scip, lambda, 0.0));

   /* for beta = - lambda / delta with delta > lambda, 
    *              f_beta   = beta - floor[beta] = 1 - ( lambda / delta ) 
    *    1 / ( 1 - f_beta )                      = delta / lambda
    */
   fbeta = 1.0 - (lambda / delta);
   onedivoneminusfbeta = delta / lambda;

   SCIPdebugMessage("     --------------------- get L1 and L2 -----------------------------------------------------\n");
   SCIPdebugMessage("     L1 = { j in N1-C1 : y*_j >= ( u_j - lambda F_{f_beta}(  u_j/delta) ) x*_j }\n");
   SCIPdebugMessage("     L2 = { j in N2-C2 : y*_j >=       - lambda F_{f_beta}(- u_j/delta)   x*_j }\n");
   
   /* set flowcover status of continuous variable x_j to 2, i.e., put j intp L1 and L2, respectively 
    *   if j is in N1\C1 and y*_j >= ( u_j - lambda F_{f_beta}(  u_j/delta) ) x*_j 
    *   if j is in N2\C2 and y*_j >=       - lambda F_{f_beta}(- u_j/delta)   x*_j 
    */
   for( j = 0; j < ntransvars; j++ ) 
   {
      SCIP_Real d;
      SCIP_Real downd;
      SCIP_Real fd;
      SCIP_Real mirval;

      assert(transvarcoefs[j] == 1 || transvarcoefs[j] == -1);
      assert(SCIPisGE(scip, transvarvubcoefs[j], 0.0));
      assert(SCIPisFeasGE(scip, transbinvarsolvals[j], 0.0) && SCIPisFeasLE(scip, transbinvarsolvals[j], 1.0));
      
      /* j in N1\C1 */
      if( transvarcoefs[j] == 1 && transvarflowcoverstatus[j] == -1 )
      {
         /* Let d = u_j/delta and alpha = f_beta, than the MIR function is defined as
          *    F_{alpha}(d) = down(d)                            , if f_d <= alpha
          *    F_{alpha}(d) = down(d) + (f_d - alpha)/(1 - alpha), if f_d >  alpha
          */ 
         d = transvarvubcoefs[j] / delta;
         downd = SCIPfloor(scip, d);
         fd = d - downd;
         if( SCIPisSumLE(scip, fd, fbeta) )
            mirval = downd; 
         else
            mirval = downd + ((fd - fbeta) * onedivoneminusfbeta);
         
         /* y*_j >= ( u_j - lambda F_{f_beta}(u_j/delta) ) x*_j */
         if( SCIPisFeasGE(scip, transcontvarsolvals[j], ( transvarvubcoefs[j] - ( lambda * mirval ) ) * transbinvarsolvals[j]) )
            transvarflowcoverstatus[j] = 2;
         
         SCIPdebugMessage("       <%d>: in N1-C1: %g ?>=?  ( %g - %g F_{f_beta}(%g)(%g) ) %g = %g  ---> fcstatus = %d\n", 
            j, transcontvarsolvals[j], transvarvubcoefs[j], lambda, d, mirval, transbinvarsolvals[j], 
            ( transvarvubcoefs[j] - ( lambda * mirval ) ) * transbinvarsolvals[j], transvarflowcoverstatus[j]);
      }

      /* j in N2\C2 */
      if( transvarcoefs[j] == -1 && transvarflowcoverstatus[j] == -1 )
      {
         /* Let d = - u_j/delta and alpha = f_beta, than the MIR function is defined as
          *    F_{alpha}(d) = down(d)                            , if f_d <= alpha
          *    F_{alpha}(d) = down(d) + (f_d - alpha)/(1 - alpha), if f_d >  alpha
          */ 
         d = - transvarvubcoefs[j] / delta;
         downd = SCIPfloor(scip, d);
         fd = d - downd;
         if( SCIPisSumLE(scip, fd, fbeta) )
            mirval = downd; 
         else
            mirval = downd + ((fd - fbeta) * onedivoneminusfbeta);
         
         /* j in N2\C2 and y*_j >= - lambda F_{f_beta}(- u_j/delta) x*_j */
         if( SCIPisFeasGE(scip, transcontvarsolvals[j], - ( lambda * mirval ) * transbinvarsolvals[j]) )
            transvarflowcoverstatus[j] = 2;
         
         SCIPdebugMessage("       <%d>: in N2-C2: %g ?>=?  - %g F_{f_beta}(-%g)(%g) %g = %g  ---> fcstatus = %d\n", 
            j, transcontvarsolvals[j], lambda, d, mirval, transbinvarsolvals[j], 
            - ( lambda * mirval ) * transbinvarsolvals[j], transvarflowcoverstatus[j]);
      }
   }
}

/** get for all problem variables with nonzero coefficient in current row the bound which should be used for the
 *  substitution routine in the c-MIR routine; this bound depends on the bound used for each variable to get associated 
 *  variable in transformed problem  
 */
static
SCIP_RETCODE getBoundsForSubstitution(
   SCIP*                 scip,               /**< SCIP data structure */ 
   SCIP_VAR**            vars,               /**< active problem variables */
   int                   nvars,              /**< number of active problem variables */
   int*                  boundsfortrans,     /**< bound used for transformation for all non-binary vars of current row */
   SCIP_BOUNDTYPE*       boundtypesfortrans, /**< type of bound used for transform. for all non-binary vars of current row */
   int*                  assoctransvars,     /**< associated var in transformed problem for all vars of current row */ 
   int*                  transvarcoefs,      /**< coefficient of all vars in transformed problem */ 
   int*                  flowcoverstatus,    /**< flow cover status of all non-binary vars in transformed problem; 
                                              *   1 if in C1 & C2, 2 if in L2, -1 N1 \ C1 & N2 \ (C2&L2) */ 
   int                   ntransvars,         /**< number of vars in transformed problem */
   int*                  boundsforsubst,     /**< pointer to store bounds that should be used for substitution in c-mir for vars */
   SCIP_BOUNDTYPE*       boundtypesforsubst  /**< pointer to store types of bounds that should be used for substitution in c-mir */
   )
{
   int j;

   assert(scip != NULL);
   assert(vars != NULL);
   assert(nvars >= 0);
   assert(boundsfortrans != NULL);
   assert(boundtypesfortrans != NULL);
   assert(assoctransvars != NULL);
   assert(transvarcoefs != NULL);
   assert(flowcoverstatus != NULL);
   assert(ntransvars >= 0 && ntransvars <= nvars);
   assert(boundsforsubst != NULL);
   assert(boundtypesforsubst != NULL);

   for( j = 0; j < nvars; j++ )
   {
      assert(SCIPvarGetProbindex(vars[j]) == j);
      assert(assoctransvars[j] == -1 || ( flowcoverstatus[assoctransvars[j]] == -1 
            || flowcoverstatus[assoctransvars[j]] == 1 || flowcoverstatus[assoctransvars[j]] == 2 ));
      
      /* variable has no associated variable in transformed problem */ 
      if( assoctransvars[j] == -1 )
      {
         assert(boundsfortrans[j] == -3);
         boundsforsubst[j] = -3;
         continue;
      }

      /* binary variable */
      if( SCIPvarIsBinary(vars[j]) )
      {
         /* j in C1 & C2 */
         if( flowcoverstatus[assoctransvars[j]] == 1 )
         {
            boundsforsubst[j] = -1;
            boundtypesforsubst[j] = SCIP_BOUNDTYPE_UPPER;
         }
         /* j in N1\C1 & N2\C2 */
         else
         {
            boundsforsubst[j] = -1;
            boundtypesforsubst[j] = SCIP_BOUNDTYPE_LOWER;
         }
      }
      /* non-binary variables */
      else
      {
         /* j in C1 & C2 & L1 & L2 */
         if( flowcoverstatus[assoctransvars[j]] >= 1 )
         {
            boundsforsubst[j] = boundsfortrans[j];
            boundtypesforsubst[j] = boundtypesfortrans[j];
         }
         /* j in N1\C1 & N2\(C2 & L2) */
         else
         {
            if( boundtypesfortrans[j] == SCIP_BOUNDTYPE_UPPER )
            {
               SCIP_Real closestlb;
               int closestlbtype;

               getClosestLb(scip, vars[j], ALLOWLOCAL, &closestlb, &closestlbtype);
               assert(!SCIPisInfinity(scip, -closestlb));
               assert(closestlbtype == -1 || closestlbtype == -2);

               boundsforsubst[j] = closestlbtype;
               boundtypesforsubst[j] = SCIP_BOUNDTYPE_LOWER;
            }      
            else
            {
               SCIP_Real closestub;
               int closestubtype;
               
               getClosestUb(scip, vars[j], ALLOWLOCAL, &closestub, &closestubtype);
               assert(!SCIPisInfinity(scip, closestub));
               assert(closestubtype == -1 || closestubtype == -2);

               boundsforsubst[j] = closestubtype;
               boundtypesforsubst[j] = SCIP_BOUNDTYPE_UPPER;
            }
         }
      }
   }

   return SCIP_OKAY;
}

/** stores nonzero elements of dense coefficient vector as sparse vector, and calculates activity and norm */
static
SCIP_RETCODE storeCutInArrays(
   SCIP*                 scip,               /**< SCIP data structure */
   int                   nvars,              /**< number of problem variables */
   SCIP_VAR**            vars,               /**< problem variables */
   SCIP_Real*            cutcoefs,           /**< dense coefficient vector */
   SCIP_Real*            varsolvals,         /**< dense variable LP solution vector */
   char                  normtype,           /**< type of norm to use for efficacy norm calculation */
   SCIP_VAR**            cutvars,            /**< array to store variables of sparse cut vector */
   SCIP_Real*            cutvals,            /**< array to store coefficients of sparse cut vector */
   int*                  cutlen,             /**< pointer to store number of nonzero entries in cut */
   SCIP_Real*            cutact,             /**< pointer to store activity of cut */
   SCIP_Real*            cutnorm             /**< pointer to store norm of cut vector */
   )
{
   SCIP_Real val;
   SCIP_Real absval;
   SCIP_Real cutsqrnorm;
   SCIP_Real act;
   SCIP_Real norm;
   int len;
   int v;

   assert(nvars == 0 || cutcoefs != NULL);
   assert(nvars == 0 || varsolvals != NULL);
   assert(cutvars != NULL);
   assert(cutvals != NULL);
   assert(cutlen != NULL);
   assert(cutact != NULL);
   assert(cutnorm != NULL);

   len = 0;
   act = 0.0;
   norm = 0.0;
   switch( normtype )
   {
   case 'e':
      cutsqrnorm = 0.0;
      for( v = 0; v < nvars; ++v )
      {
         val = cutcoefs[v];
         if( !SCIPisZero(scip, val) )
         {
            act += val * varsolvals[v];
            cutsqrnorm += SQR(val);
            cutvars[len] = vars[v];
            cutvals[len] = val;
            len++;
         }
      }
      norm = SQRT(cutsqrnorm);
      break;
   case 'm':
      for( v = 0; v < nvars; ++v )
      {
         val = cutcoefs[v];
         if( !SCIPisZero(scip, val) )
         {
            act += val * varsolvals[v];
            absval = REALABS(val);
            norm = MAX(norm, absval);
            cutvars[len] = vars[v];
            cutvals[len] = val;
            len++;
         }
      }
      break;
   case 's':
      for( v = 0; v < nvars; ++v )
      {
         val = cutcoefs[v];
         if( !SCIPisZero(scip, val) )
         {
            act += val * varsolvals[v];
            norm += REALABS(val);
            cutvars[len] = vars[v];
            cutvals[len] = val;
            len++;
         }
      }
      break;
   case 'd':
      for( v = 0; v < nvars; ++v )
      {
         val = cutcoefs[v];
         if( !SCIPisZero(scip, val) )
         {
            act += val * varsolvals[v];
            norm = 1.0;
            cutvars[len] = vars[v];
            cutvals[len] = val;
            len++;
         }
      }
      break;
   default:
      SCIPerrorMessage("invalid efficacy norm parameter '%c'\n", normtype);
      return SCIP_INVALIDDATA;
   }

   *cutlen = len;
   *cutact = act;
   *cutnorm = norm;

   return SCIP_OKAY;
}

/** adds given cut to LP if violated */
static
SCIP_RETCODE addCut(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SEPA*            sepa,               /**< separator */
   SCIP_SEPADATA*        sepadata,           /**< separator data */
   SCIP_VAR**            vars,               /**< problem variables */
   int                   nvars,              /**< number of problem variables */
   SCIP_SOL*             sol,                /**< the solution that should be separated, or NULL for LP solution */
   SCIP_Real*            varsolvals,         /**< solution values of active variables */
   SCIP_Real*            cutcoefs,           /**< coefficients of active variables in cut */
   SCIP_Real             cutrhs,             /**< right hand side of cut */
   SCIP_Bool             cutislocal,         /**< is the cut only locally valid? */
   char                  normtype,           /**< type of norm to use for efficacy norm calculation */
   int*                  ncuts               /**< pointer to count the number of added cuts */
   )
{
   SCIP_VAR** cutvars;
   SCIP_Real* cutvals;
   SCIP_Real cutact;
   SCIP_Real cutnorm;
   int cutlen;
   SCIP_Bool success;

   assert(scip != NULL);
   assert(varsolvals != NULL);
   assert(cutcoefs != NULL);
   assert(ncuts != NULL);
   assert(nvars == 0 || vars != NULL);

   SCIPdebugMessage("--------------------- found cut ---------------------------------------------------------\n");

   /* gets temporary memory for storing the cut as sparse row */
   SCIP_CALL( SCIPallocBufferArray(scip, &cutvars, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &cutvals, nvars) );

   /* stores the cut as sparse row, calculates activity and norm of cut */
   SCIP_CALL( storeCutInArrays(scip, nvars, vars, cutcoefs, varsolvals, normtype,
         cutvars, cutvals, &cutlen, &cutact, &cutnorm) );

   if( SCIPisPositive(scip, cutnorm) && SCIPisEfficacious(scip, (cutact - cutrhs)/cutnorm) )
   {
      SCIP_ROW* cut;
      char cutname[SCIP_MAXSTRLEN];

      /* creates the cut */
      (void) SCIPsnprintf(cutname, SCIP_MAXSTRLEN, "flowcover%d_%d", SCIPgetNLPs(scip), *ncuts);
      SCIP_CALL( SCIPcreateEmptyRowSepa(scip, &cut, sepa, cutname, -SCIPinfinity(scip), cutrhs,
            cutislocal, FALSE, sepadata->dynamiccuts) );
      SCIP_CALL( SCIPaddVarsToRow(scip, cut, cutlen, cutvars, cutvals) );

      SCIPdebugMessage(" -> found potential flowcover cut <%s>: activity=%f, rhs=%f, norm=%f, eff=%f\n",
         cutname, cutact, cutrhs, cutnorm, SCIPgetCutEfficacy(scip, sol, cut));
      SCIPdebug( SCIP_CALL( SCIPprintRow(scip, cut, NULL) ) );

#if 0 /* tries to scale the cut to integral values */
      SCIP_CALL( SCIPmakeRowIntegral(scip, cut, -SCIPepsilon(scip), SCIPsumepsilon(scip),
            10, 100.0, MAKECONTINTEGRAL, &success) );
      if( success && !SCIPisCutEfficacious(scip, sol, cut) )
      {
         SCIPdebugMessage(" -> flowcover cut <%s> no longer efficacious: act=%f, rhs=%f, norm=%f, eff=%f\n",
            cutname, cutact, cutrhs, cutnorm, SCIPgetCutEfficacy(scip, sol, cut));
         SCIPdebug( SCIP_CALL( SCIPprintRow(scip, cut, NULL) ) );
         success = FALSE;
      }
#else
      success = TRUE;
#endif

      /* if scaling was successful, adds the cut */
      if( success ) /*lint !e774*/ /* Boolean within 'if' always evaluates to True */
      {
         SCIPdebugMessage(" -> found flowcover cut <%s>: act=%f, rhs=%f, norm=%f, eff=%f, min=%f, max=%f (range=%g)\n",
            cutname, cutact, cutrhs, cutnorm, SCIPgetCutEfficacy(scip, sol, cut),
            SCIPgetRowMinCoef(scip, cut), SCIPgetRowMaxCoef(scip, cut),
            SCIPgetRowMaxCoef(scip, cut)/SCIPgetRowMinCoef(scip, cut));
         SCIPdebug( SCIP_CALL( SCIPprintRow(scip, cut, NULL) ) );
         SCIP_CALL( SCIPaddCut(scip, sol, cut, FALSE) );
         if( !cutislocal )
         {
            SCIP_CALL( SCIPaddPoolCut(scip, cut) );
         }
         (*ncuts)++;
      }
      
      /* releases the row */
      SCIP_CALL( SCIPreleaseRow(scip, &cut) );
   }
   
   /* frees temporary memory */
   SCIPfreeBufferArray(scip, &cutvals);
   SCIPfreeBufferArray(scip, &cutvars);

   return SCIP_OKAY;   
}

/** calculates efficacy of the given cut */
static
SCIP_Real calcEfficacy(
   int                   nvars,              /**< number of variables in the problem */
   SCIP_Real*            cutcoefs,           /**< dense vector of cut coefficients */
   SCIP_Real             cutrhs,             /**< right hand side of cut */
   SCIP_Real             cutact              /**< activity of cut */
   )
{
   SCIP_Real sqrnorm;
   int i;
   
   assert(cutcoefs != NULL);

   sqrnorm = 0;
   for( i = 0; i < nvars; ++i )
      sqrnorm += SQR(cutcoefs[i]);
   sqrnorm = MAX(sqrnorm, 1e-06);
   return (cutact - cutrhs)/SQRT(sqrnorm);
}

/* generate c-MIRFCI for different sets L1 and L2 and different values of delta */
static
SCIP_RETCODE cutGenerationHeuristic(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SEPA*            sepa,               /**< separator */
   SCIP_SEPADATA*        sepadata,           /**< separator data */
   SCIP_VAR**            vars,               /**< active problem variables */
   int                   nvars,              /**< number of active problem variables */
   SCIP_SOL*             sol,                /**< the solution that should be separated, or NULL for LP solution */
   SCIP_Real*            varsolvals,         /**< solution values of active variables */
   SCIP_Real*            rowweights,         /**< weight of rows in aggregated row */ 
   SCIP_Real             scalar,             /**< additional scaling factor of rows in aggregation */
   int*                  boundsfortrans,     /**< bound used for all non-bin vars of row */
   SCIP_BOUNDTYPE*       boundtypesfortrans, /**< type of bound used for all non-binary vars of row */
   int*                  assoctransvars,     /**< associated var in relaxed set for all vars of row */ 
   int                   ntransvars,         /**< number of real variables in N1&N2 */
   int*                  transvarcoefs,      /**< coefficient of all continuous variables in N1 & N2 */ 
   SCIP_Real*            transbinvarsolvals, /**< LP solution value of binary variable in vub of all real vars in N1&N2 */
   SCIP_Real*            transcontvarsolvals,/**< LP solution value of all real vars in N1&N2 */
   SCIP_Real*            transvarvubcoefs,   /**< coefficient of vub of all continuous variables in N1 & N2 */
   int*                  transvarflowcoverstatus, /**< pointer to store whether non-binary var is in L2 (2) or not (-1 or 1) */ 
   SCIP_Real             lambda,             /**< lambda */
   char                  normtype,           /**< type of norm to use for efficacy norm calculation */
   int*                  ncuts               /**< pointer to count the number of generated cuts */
   )
{
   SCIP_BOUNDTYPE* boundtypesforsubst;
   SCIP_Real* candsetdelta;
   SCIP_Real* cutcoefs;      
   SCIP_Real* testeddeltas;
   int* boundsforsubst;
   int* transvarflowcoverstatustmp;
   SCIP_Real bestdelta;
   SCIP_Real bestefficacy;
   SCIP_Real c1vubcoefsmax;
   SCIP_Real cutrhs;         
   SCIP_Real cutact;
   SCIP_Real l2tmpvubcoefsmax;
   SCIP_Real nvubcoefsmax;
   SCIP_Bool cutislocal;
   SCIP_Bool success;
   int ncandsetdelta;
   int ntesteddeltas;
   int startidx;
   int j;

   assert(scip != NULL);
   assert(sepadata != NULL);
   assert(vars != NULL);
   assert(varsolvals != NULL);
   assert(rowweights != NULL);
   assert(boundsfortrans != NULL);
   assert(boundtypesfortrans != NULL);
   assert(assoctransvars != NULL);
   assert(ntransvars >=0);
   assert(transvarcoefs != NULL);
   assert(transbinvarsolvals != NULL);
   assert(transcontvarsolvals != NULL);
   assert(transvarvubcoefs != NULL);
   assert(transvarflowcoverstatus != NULL);
   assert(SCIPisFeasGT(scip, lambda, 0.0));
   assert(ncuts != NULL);
      
   SCIPdebugMessage("--------------------- cut generation heuristic ------------------------------------------\n");

   /* get data structures */
   SCIP_CALL( SCIPallocBufferArray(scip, &testeddeltas, ntransvars + 4) );
   SCIP_CALL( SCIPallocBufferArray(scip, &candsetdelta, ntransvars + 4) );
   SCIP_CALL( SCIPallocBufferArray(scip, &transvarflowcoverstatustmp, ntransvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &boundsforsubst, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &boundtypesforsubst, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &cutcoefs, nvars) );

   SCIPdebugMessage("1. get candidate set for the value of delta:\n");

   /* get candidate set for the value of delta  
    *   N* = { u_j : j in N and u_j > lambda} & { max{ u_j : j in N and u_j >= lambda } + 1, lambda + 1 };
    * store the following values which will be tested in any case at the beginning of the array
    *   max{ u_j : j in N   and u_j >= lambda } + 1, 
    *                                  lambda   + 1, 
    *   max{ u_j : j in C1  and u_j >  lambda }, 
    *   max{ u_j : j in L~2 and u_j >  lambda } 
    */
   ncandsetdelta = 0;
   startidx = 4;

   /* get max{ u_j : j in C1  and u_j > lambda }, max{ u_j : j in L~2 and u_j > lambda }, and 
    * max{ u_j : j in N and u_j >= lambda } and store { u_j : j in N and u_j > lambda } 
    */
   c1vubcoefsmax = SCIP_REAL_MIN;
   l2tmpvubcoefsmax = SCIP_REAL_MIN;
   nvubcoefsmax = SCIP_REAL_MIN;
   for( j = 0; j < ntransvars; j++ )
   {
      SCIP_Real val;

      /* j is in C1 and u_j > lambda */
      if( transvarcoefs[j] == 1 && transvarflowcoverstatus[j] == 1 
         && SCIPisFeasGT(scip, transvarvubcoefs[j], lambda) )
         c1vubcoefsmax = MAX(c1vubcoefsmax, transvarvubcoefs[j]);
         
      /* j is in L~2 and u_j > lambda */
      val = MIN(transvarvubcoefs[j], lambda) * transbinvarsolvals[j];
      if( transvarcoefs[j] == -1 && transvarflowcoverstatus[j] == -1 && SCIPisFeasLE(scip, val, transcontvarsolvals[j])
         && SCIPisFeasGT(scip, transvarvubcoefs[j], lambda) )
         l2tmpvubcoefsmax = MAX(l2tmpvubcoefsmax, transvarvubcoefs[j]);

      /* u_j + 1 > lambda */
      if( SCIPisFeasGT(scip, transvarvubcoefs[j] + 1.0, lambda) )
         nvubcoefsmax = MAX(nvubcoefsmax, transvarvubcoefs[j]);
         
      /* u_j > lambda */
      if( SCIPisFeasGT(scip, transvarvubcoefs[j], lambda) )
      {
         candsetdelta[startidx + ncandsetdelta] = transvarvubcoefs[j];
         ncandsetdelta++;
         SCIPdebugMessage("     u_j              = %g\n", transvarvubcoefs[j]);
      }
   }
   
   /* store max { u_j : j in N and u_j >= lambda } + 1 */
   if( !SCIPisInfinity(scip, -nvubcoefsmax) )
   {
      assert(SCIPisFeasGT(scip, nvubcoefsmax + 1.0, lambda));
      startidx--;
      candsetdelta[startidx] = nvubcoefsmax + 1.0;
      ncandsetdelta++;
      SCIPdebugMessage("     max u_j+1        = %g\n", nvubcoefsmax + 1.0);
   }

   /* store lambda + 1 */
   startidx--;
   candsetdelta[startidx] = lambda + 1.0;
   ncandsetdelta++;
   SCIPdebugMessage("     lambda + 1       = %g\n", lambda + 1.0);

   /* store max{ u_j : j in C1 and u_j > lambda } and max{ u_j : j in C1 & L~2 and u_j > lambda } */
   if( !SCIPisInfinity(scip, -l2tmpvubcoefsmax) && SCIPisFeasGT(scip, l2tmpvubcoefsmax, c1vubcoefsmax) )
   {
      assert(SCIPisFeasGT(scip, l2tmpvubcoefsmax, lambda));
      startidx--;
      candsetdelta[startidx] = l2tmpvubcoefsmax;
      ncandsetdelta++;
      SCIPdebugMessage("     l2tmpvubcoefsmax = %g\n", l2tmpvubcoefsmax);
   }
   if( !SCIPisInfinity(scip, -c1vubcoefsmax) )
   {
      assert(SCIPisFeasGT(scip, c1vubcoefsmax, lambda));
      startidx--;
      candsetdelta[startidx] = c1vubcoefsmax;
      ncandsetdelta++;
      SCIPdebugMessage("     c1vubcoefsmax    = %g\n", c1vubcoefsmax);
   }
     
   assert(startidx >= 0 && startidx <= 3);
   assert(startidx + ncandsetdelta >= 4);
   assert(ncandsetdelta >= 1 && ncandsetdelta <= ntransvars + 4);

   SCIPdebugMessage("2. generate c-MIRFCIs for different values of delta:\n");

   /* for each value of delta choose L1 subset N1\C1 and L2 subset N2\C2 by comparison, generate the 
    * c-MIRFCI for delta, (C1, C2) and (L1, L2) and select the most efficient c-MIRFCI
    */ 
   ntesteddeltas = 0;
   bestdelta = 0.0;
   bestefficacy = 0.0;
   for( j = startidx; j < startidx + ncandsetdelta && ntesteddeltas < sepadata->maxtestdelta; j++ )
   {
      SCIP_Real delta;
      SCIP_Real onedivdelta;
      SCIP_Bool tested;
      int i;

      /* get current delta and corresponding scalar for c-MIR routine */
      delta = candsetdelta[j];
      onedivdelta = 1.0 / delta;

      /* do not use scaling factors (1/delta) which are too small, since this max cause numerical problems;
       * besides for relatively small lambda and large scaling factors (delta), we get f_beta = 1 - lambda/delta > MINFRAC 
       */
      if( SCIPisZero(scip, scalar * onedivdelta) )
         continue;
         
      /* check, if delta was already tested */
      tested = FALSE;
      for( i = 0; i < ntesteddeltas && !tested; i++ )
         tested = SCIPisEQ(scip, testeddeltas[i], delta);
      if( tested )
         continue;
      testeddeltas[ntesteddeltas] = delta;
      ntesteddeltas++;

      SCIPdebugMessage("   delta = %g:\n", delta);

      /* work on copy of transvarflowcoverstatus because current choice of sets L1 and L2 will change 
       * transvarflowcoverstatus 
       */  
      BMScopyMemoryArray(transvarflowcoverstatustmp, transvarflowcoverstatus, ntransvars);

      /* get L1 subset of N1\C1 and L2 subset of N2\C2 by comparison */
      getL1L2(scip, ntransvars, transvarcoefs, transbinvarsolvals, transcontvarsolvals, transvarvubcoefs, 
         transvarflowcoverstatustmp, delta, lambda);
      
      /* get bounds for substitution in c-MIR routine for original mixed integer set;
       * note that the scalar has already been considered in the constructed 0-1 single node flow relaxation 
       */
      SCIP_CALL( getBoundsForSubstitution(scip, vars, nvars, boundsfortrans, boundtypesfortrans, assoctransvars, 
            transvarcoefs, transvarflowcoverstatustmp, ntransvars, boundsforsubst, boundtypesforsubst ) );
      
      /* generate c-MIRFCI for flow cover (C1,C2), L1 subset N1\C1 and L2 subset N2\C2 and delta */
      SCIP_CALL( SCIPcalcMIR(scip, sol, BOUNDSWITCH, TRUE, ALLOWLOCAL, FIXINTEGRALRHS, boundsforsubst, boundtypesforsubst,
            (int) MAXAGGRLEN(nvars), 1.0, MINFRAC, MAXFRAC, rowweights, scalar * onedivdelta, NULL, NULL, cutcoefs, 
            &cutrhs, &cutact, &success, &cutislocal) );
      assert(ALLOWLOCAL || !cutislocal);
      
      /* delta leads to c-MIRFCI which is more violated */
      if( success )
      {
         SCIP_Real efficacy;
         
         for(i = 0; i < nvars; i++)
            cutcoefs[i] = lambda * cutcoefs[i];
         cutrhs = lambda * cutrhs;
         cutact = lambda * cutact;
         
         efficacy = calcEfficacy(nvars, cutcoefs, cutrhs, cutact);
         
         SCIPdebugMessage("   ---> act = %g  rhs = %g  eff = %g (old besteff = %g, old bestdelta=%g)\n", 
            cutact, cutrhs, efficacy, bestefficacy, bestdelta);

         if( efficacy > bestefficacy )
         {
            bestdelta = delta;
            bestefficacy = efficacy;
         }
      }
   }

   /* delta found */
   if( SCIPisEfficacious(scip, bestefficacy) )
   {
      SCIP_Real onedivbestdelta;
      int i;

      assert(bestdelta != 0.0);
      onedivbestdelta = 1.0 / bestdelta;

      /* for best value of delta: get L1 subset of N1\C1 and L2 subset of N2\C2 by comparison */
      getL1L2(scip, ntransvars, transvarcoefs, transbinvarsolvals, transcontvarsolvals, transvarvubcoefs, 
         transvarflowcoverstatus, bestdelta, lambda);
      
      /* for best value of delta: get bounds for substitution in c-MIR routine for original mixed integer set 
       * note that the scalar has already been considered in the constructed 0-1 single node flow relaxation 
       */
      SCIP_CALL( getBoundsForSubstitution(scip, vars, nvars, boundsfortrans, boundtypesfortrans, assoctransvars, 
            transvarcoefs, transvarflowcoverstatus, ntransvars, boundsforsubst, boundtypesforsubst ) );
      
      /* generate c-MIRFCI for flow cover (C1,C2), L1 subset N1\C1 and L2 subset N2\C2 and bestdelta */
      SCIP_CALL( SCIPcalcMIR(scip, sol, BOUNDSWITCH, TRUE, ALLOWLOCAL, FIXINTEGRALRHS, boundsforsubst, boundtypesforsubst,
            (int) MAXAGGRLEN(nvars), 1.0, MINFRAC, MAXFRAC, rowweights, scalar * onedivbestdelta, NULL, NULL, cutcoefs, 
            &cutrhs, &cutact, &success, &cutislocal) );
      assert(ALLOWLOCAL || !cutislocal);
      assert(success); 
      
      for(i = 0; i < nvars; i++)
         cutcoefs[i] = lambda * cutcoefs[i];
      cutrhs = lambda * cutrhs;
      cutact = lambda * cutact;

      assert(SCIPisFeasEQ(scip, bestefficacy, calcEfficacy(nvars, cutcoefs, cutrhs, cutact)));
      SCIP_CALL( addCut(scip, sepa, sepadata, vars, nvars, sol, varsolvals, cutcoefs, cutrhs, cutislocal, normtype, ncuts) );
   }

   /* free data structures */
   SCIPfreeBufferArray(scip, &cutcoefs);
   SCIPfreeBufferArray(scip, &boundtypesforsubst);
   SCIPfreeBufferArray(scip, &boundsforsubst);
   SCIPfreeBufferArray(scip, &transvarflowcoverstatustmp);
   SCIPfreeBufferArray(scip, &candsetdelta);
   SCIPfreeBufferArray(scip, &testeddeltas);
   
   return SCIP_OKAY;
}

/** search and add flowcover cuts that separate the given primal solution */
static
SCIP_RETCODE separateCuts(
   SCIP*                 scip,               /**< SCIP data structure */ 
   SCIP_SEPA*            sepa,               /**< the flowcover separator */
   SCIP_SOL*             sol,                /**< the solution that should be separated, or NULL for LP solution */
   SCIP_RESULT*          result              /**< pointer to store the result */
   )
{
   SCIP_ROW** rows;     
   SCIP_VAR** vars;
   SCIP_BOUNDTYPE* boundtypesfortrans;
   SCIP_SEPADATA* sepadata;
   SCIP_Real* rowlhsscores;
   SCIP_Real* rowrhsscores;
   SCIP_Real* rowscores;
   SCIP_Real* rowweights;  
   SCIP_Real* transbinvarsolvals;
   SCIP_Real* transcontvarsolvals;
   SCIP_Real* transvarvubcoefs;
   SCIP_Real* varsolvals;
   int* assoctransvars;
   int* boundsfortrans;
   int* covervars;
   int* noncovervars;
   int* roworder;
   int* transvarcoefs;
   int* transvarflowcoverstatus;
   SCIP_Real lambda;
   SCIP_Real maxslack;
   SCIP_Real objnorm;
   SCIP_Real transcapacity;
   SCIP_Bool transsuccess;
   SCIP_Bool flowcoverfound;
   char normtype;
   int depth;
   int maxfails;
   int maxsepacuts;
   int maxtries;
   int ncalls;
   int ncovervars;
   int nfails;
   int nnoncovervars;
   int nrows;
   int ntransvars;
   int ntries;
   int nvars;
   int ncuts;
   int r;

   assert(scip != NULL);
   assert(sepa != NULL);
   assert(result != NULL);
   assert(*result == SCIP_DIDNOTRUN);
 
   sepadata = SCIPsepaGetData(sepa);
   assert(sepadata != NULL);

   depth = SCIPgetDepth(scip);
   ncalls = SCIPsepaGetNCallsAtNode(sepa);

   /* only call the flow cover cuts separator a given number of times at each node */
   if( (depth == 0 && sepadata->maxroundsroot >= 0 && ncalls >= sepadata->maxroundsroot)
      || (depth > 0 && sepadata->maxrounds >= 0 && ncalls >= sepadata->maxrounds) )
      return SCIP_OKAY;

   /* get all rows and number of rows */
   SCIP_CALL( SCIPgetLPRowsData(scip, &rows, &nrows) ); 
   assert(nrows == 0 || rows != NULL);

   /* nothing to do, if LP is empty */
   if( nrows == 0 )
      return SCIP_OKAY;

   /* get active problem variables */
   vars = SCIPgetVars(scip);
   nvars = SCIPgetNVars(scip);
   assert(nvars == 0 || vars != NULL);

   /* nothing to do, if problem has no variables */
   if( nvars == 0 )
      return SCIP_OKAY;

   /* check whether SCIP was stopped in the meantime */
   if( SCIPisStopped(scip) )
      return SCIP_OKAY;
   
   *result = SCIP_DIDNOTFIND;

   /* get the type of norm to use for efficacy calculations */
   SCIP_CALL( SCIPgetCharParam(scip, "separating/efficacynorm", &normtype) );

   /* get data structures */
   SCIP_CALL( SCIPallocBufferArray(scip, &rowlhsscores, nrows) );
   SCIP_CALL( SCIPallocBufferArray(scip, &rowrhsscores, nrows) );
   SCIP_CALL( SCIPallocBufferArray(scip, &rowscores, nrows) );
   SCIP_CALL( SCIPallocBufferArray(scip, &roworder, nrows) );
   SCIP_CALL( SCIPallocBufferArray(scip, &varsolvals, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &assoctransvars, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &transvarcoefs, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &transvarvubcoefs, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &transbinvarsolvals, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &transcontvarsolvals, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &transvarflowcoverstatus, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &boundsfortrans, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &boundtypesfortrans, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &covervars, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &noncovervars, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &rowweights, nrows) );

   /* get the solution values for all active variables */
   SCIP_CALL( SCIPgetSolVals(scip, sol, nvars, vars, varsolvals) );
   
   /* get the maximal number of cuts allowed in a separation round */
   if( depth == 0 )
   {
      maxtries = sepadata->maxtriesroot;
      maxfails = sepadata->maxfailsroot;
      maxsepacuts = sepadata->maxsepacutsroot;
      maxslack = sepadata->maxslackroot;
   }   
   else
   {
      maxtries = sepadata->maxtries;
      maxfails = sepadata->maxfails;
      maxsepacuts = sepadata->maxsepacuts;
      maxslack = sepadata->maxslack;
   }

   /* calculate row scores for both sides of all rows, and sort rows by nonincreasing maximal score */
   objnorm = SCIPgetObjNorm(scip);
   objnorm = MAX(objnorm, 1.0);
   for( r = 0; r < nrows; r++ )
   {
      int nnonz;
      int i;

      assert(SCIProwGetLPPos(rows[r]) == r);

      nnonz = SCIProwGetNLPNonz(rows[r]);
      if( nnonz == 0 )
      {
         /* ignore empty rows */
         rowlhsscores[r] = 0.0;
         rowrhsscores[r] = 0.0;
      }
      else
      {
         SCIP_Real activity;
         SCIP_Real lhs;
         SCIP_Real rhs;
         SCIP_Real dualsol;
         SCIP_Real dualscore;
         SCIP_Real rowdensity;
         SCIP_Real rownorm;
         SCIP_Real slack;

         dualsol = (sol == NULL ? SCIProwGetDualsol(rows[r]) : 1.0);
         activity = SCIPgetRowSolActivity(scip, rows[r], sol);
         lhs = SCIProwGetLhs(rows[r]);
         rhs = SCIProwGetRhs(rows[r]);
         rownorm = SCIProwGetNorm(rows[r]);
         rownorm = MAX(rownorm, 0.1);
         rowdensity = (SCIP_Real)nnonz/(SCIP_Real)nvars;
         assert(SCIPisPositive(scip, rownorm));

         slack = (activity - lhs)/rownorm;
         dualscore = MAX(dualsol/objnorm, 0.0001);
         if( !SCIPisInfinity(scip, -lhs) && SCIPisLE(scip, slack, maxslack)
            && (ALLOWLOCAL || !SCIProwIsLocal(rows[r])) /*lint !e506 !e774*/
            && rowdensity <= sepadata->maxrowdensity )  
         {
            rowlhsscores[r] = dualscore + DENSSCORE * rowdensity + sepadata->slackscore * MAX(1.0 - slack, 0.0);
            assert(rowlhsscores[r] > 0.0);
         }
         else
            rowlhsscores[r] = 0.0;

         slack = (rhs - activity)/rownorm;
         dualscore = MAX(-dualsol/objnorm, 0.0001);
         if( !SCIPisInfinity(scip, rhs) && SCIPisLE(scip, slack, maxslack)
            && (ALLOWLOCAL || !SCIProwIsLocal(rows[r])) /*lint !e506 !e774*/
            && rowdensity <= sepadata->maxrowdensity )
         {
            rowrhsscores[r] = dualscore + DENSSCORE * rowdensity + sepadata->slackscore * MAX(1.0 - slack, 0.0); 
            assert(rowrhsscores[r] > 0.0);
         }
         else
            rowrhsscores[r] = 0.0;
      }
      rowscores[r] = MAX(rowlhsscores[r], rowrhsscores[r]);
      for( i = r; i > 0 && rowscores[r] > rowscores[roworder[i-1]]; --i )
         roworder[i] = roworder[i-1];
      assert(0 <= i && i <= r);
      roworder[i] = r;
   }

   /* initialize weights of rows for aggregation for c-mir routine */
   BMSclearMemoryArray(rowweights, nrows);

   /* try to generate a flow cover cut for each row in the LP */
   ncuts = 0;
   if( maxtries < 0 )
      maxtries = INT_MAX;
   if( maxfails < 0 )
      maxfails = INT_MAX;
   else if( depth == 0 && 2*SCIPgetNSepaRounds(scip) < maxfails )
      maxfails += maxfails - 2*SCIPgetNSepaRounds(scip); /* allow up to double as many fails in early separounds of root node */
   ntries = 0;
   nfails = 0;
   for( r = 0; r < nrows && ntries < maxtries && ncuts < maxsepacuts && rowscores[roworder[r]] > 0.0
           && !SCIPisStopped(scip); r++ )
   {
      SCIP_Bool wastried;
      int oldncuts;
      SCIP_Real rowact;
      SCIP_Real mult;
      
      oldncuts = ncuts;
      wastried = FALSE;

      /* update weight of rows for aggregation in c-MIR routine; all rows but current one have weight 0.0 */
      if( r > 0 )
      {
         assert(rowweights[roworder[r-1]] == 1.0 || rowweights[roworder[r-1]] == -1.0);
         rowweights[roworder[r-1]] = 0.0;
      }

      /* decide which side of the row should be used */
      rowact = SCIPgetRowSolActivity(scip, rows[roworder[r]], sol);
      if( rowact < 0.5 * SCIProwGetLhs(rows[roworder[r]]) + 0.5 * SCIProwGetRhs(rows[roworder[r]]) )
         rowweights[roworder[r]] = -1.0;
      else 
         rowweights[roworder[r]] = 1.0;

      SCIPdebugMessage("===================== flow cover separation for row <%s> (%d of %d) ===================== \n",
         SCIProwGetName(rows[roworder[r]]), r, nrows);
      SCIPdebug( SCIP_CALL( SCIPprintRow(scip, rows[roworder[r]], NULL) ) );
      SCIPdebugMessage("rowact=%g is closer to %s --> rowweight=%g\n", rowact, 
         rowweights[roworder[r]] == 1 ? "rhs" : "lhs", rowweights[roworder[r]]);

      mult = 1.0;
      do
      {
         assert(mult == 1.0 || (mult == -1.0 && sepadata->multbyminusone));

         /* construct 0-1 single node flow relaxation (with some additional simple constraints) of the mixed integer set 
          * corresponding to the current row  
          *       sum_{j in N} a_j x_j + sum_{j in M} c_j y_j   + s =     rhs + const   or
          *   - ( sum_{j in N} a_j x_j + sum_{j in M} c_j y_j ) + s = - ( lhs - const )
          * multiplied by mult in { +1, -1 }
          */
         SCIP_CALL( constructSNFRelaxation(scip, vars, nvars, varsolvals, rows[roworder[r]], rowweights[roworder[r]], 
               mult, boundsfortrans, boundtypesfortrans, assoctransvars, transvarcoefs, transbinvarsolvals, 
               transcontvarsolvals, transvarvubcoefs, &ntransvars, &transcapacity, &transsuccess) );
         if( !transsuccess )
         {
            /* transformation will fail for mult = -1, too */
            SCIPdebugMessage("mult=%g: no 0-1 single node flow relaxation found\n", mult);
            break;
         }  
       
	 flowcoverfound = FALSE;

         /* get a flow cover (C1, C2) for the constructed 0-1 single node flow set */
         SCIP_CALL( getFlowCover(scip, transvarcoefs, transbinvarsolvals, transvarvubcoefs, ntransvars, transcapacity, 
               &ncovervars, &nnoncovervars, transvarflowcoverstatus, &lambda, &flowcoverfound) );
         if( !flowcoverfound ) 
         {
            SCIPdebugMessage("mult=%g: no flow cover found\n", mult);
            mult *= -1.0;
            continue; 
         }
         assert(SCIPisFeasGT(scip, lambda, 0.0)); 

         /* generate most violated c-MIRFCI for different sets L1 and L2 and different values of delta and add it to the LP */
         SCIP_CALL( cutGenerationHeuristic(scip, sepa, sepadata, vars, nvars, sol, varsolvals, rowweights, mult, boundsfortrans,
               boundtypesfortrans, assoctransvars, ntransvars, transvarcoefs, transbinvarsolvals, transcontvarsolvals,
               transvarvubcoefs, transvarflowcoverstatus, lambda, normtype, &ncuts) );

         wastried = TRUE;
         mult *= -1.0;
      }
      while( sepadata->multbyminusone && mult < 0.0 );

      if( !wastried )
         continue;

      ntries++;
      if( ncuts == oldncuts )
      {
         nfails++;
         if( nfails >= maxfails )
            break;
      }
      else
         nfails = 0;
   }   

   /* free data structures */
   SCIPfreeBufferArray(scip, &rowweights);
   SCIPfreeBufferArray(scip, &noncovervars);
   SCIPfreeBufferArray(scip, &covervars);
   SCIPfreeBufferArray(scip, &boundtypesfortrans);
   SCIPfreeBufferArray(scip, &boundsfortrans);
   SCIPfreeBufferArray(scip, &transvarflowcoverstatus);
   SCIPfreeBufferArray(scip, &transcontvarsolvals);
   SCIPfreeBufferArray(scip, &transbinvarsolvals);
   SCIPfreeBufferArray(scip, &transvarvubcoefs);
   SCIPfreeBufferArray(scip, &transvarcoefs);
   SCIPfreeBufferArray(scip, &assoctransvars);
   SCIPfreeBufferArray(scip, &varsolvals);
   SCIPfreeBufferArray(scip, &roworder);
   SCIPfreeBufferArray(scip, &rowscores);
   SCIPfreeBufferArray(scip, &rowrhsscores);
   SCIPfreeBufferArray(scip, &rowlhsscores);

   if( ncuts > 0 )
      *result = SCIP_SEPARATED;
    
   return SCIP_OKAY;
}


/*
 * Callback methods of separator
 */

/** copy method for separator plugins (called when SCIP copies plugins) */
static
SCIP_DECL_SEPACOPY(sepaCopyFlowcover)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(sepa != NULL);
   assert(strcmp(SCIPsepaGetName(sepa), SEPA_NAME) == 0);

   /* call inclusion method of constraint handler */
   SCIP_CALL( SCIPincludeSepaFlowcover(scip) );
 
   return SCIP_OKAY;
}

/** destructor of separator to free user data (called when SCIP is exiting) */
static
SCIP_DECL_SEPAFREE(sepaFreeFlowcover)
{  /*lint --e{715}*/
   SCIP_SEPADATA* sepadata;

   /* free separator data */
   sepadata = SCIPsepaGetData(sepa);
   assert(sepadata != NULL);

   SCIPfreeMemory(scip, &sepadata);

   SCIPsepaSetData(sepa, NULL);

   return SCIP_OKAY;
}


/** LP solution separation method of separator */
static
SCIP_DECL_SEPAEXECLP(sepaExeclpFlowcover)
{  /*lint --e{715}*/

   *result = SCIP_DIDNOTRUN;

   /* only call separator, if we are not close to terminating */
   if( SCIPisStopped(scip) )
      return SCIP_OKAY;

   /* only call separator, if an optimal LP solution is at hand */
   if( SCIPgetLPSolstat(scip) != SCIP_LPSOLSTAT_OPTIMAL )
      return SCIP_OKAY;

   /* only call separator, if there are fractional variables */
   if( SCIPgetNLPBranchCands(scip) == 0 )
      return SCIP_OKAY;

   SCIP_CALL( separateCuts(scip, sepa, NULL, result) );

   return SCIP_OKAY;
}


/** arbitrary primal solution separation method of separator */
static
SCIP_DECL_SEPAEXECSOL(sepaExecsolFlowcover)
{  /*lint --e{715}*/

   *result = SCIP_DIDNOTRUN;

   SCIP_CALL( separateCuts(scip, sepa, sol, result) );

   return SCIP_OKAY;
}


/*
 * separator specific interface methods
 */

/** creates the flowcover separator and includes it in SCIP */
SCIP_RETCODE SCIPincludeSepaFlowcover(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_SEPADATA* sepadata;
   SCIP_SEPA* sepa;

   /* create flowcover separator data */
   SCIP_CALL( SCIPallocMemory(scip, &sepadata) );

   /* include separator */
   SCIP_CALL( SCIPincludeSepaBasic(scip, &sepa, SEPA_NAME, SEPA_DESC, SEPA_PRIORITY, SEPA_FREQ, SEPA_MAXBOUNDDIST,
         SEPA_USESSUBSCIP, SEPA_DELAY,
         sepaExeclpFlowcover, sepaExecsolFlowcover,
         sepadata) );

   assert(sepa != NULL);

   /* set non-NULL pointers to callback methods */
   SCIP_CALL( SCIPsetSepaCopy(scip, sepa, sepaCopyFlowcover) );
   SCIP_CALL( SCIPsetSepaFree(scip, sepa, sepaFreeFlowcover) );

   /* add flow cover cuts separator parameters */
   SCIP_CALL( SCIPaddIntParam(scip,
         "separating/flowcover/maxrounds",
         "maximal number of separation rounds per node (-1: unlimited)",
         &sepadata->maxrounds, FALSE, DEFAULT_MAXROUNDS, -1, INT_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddIntParam(scip,
         "separating/flowcover/maxroundsroot",
         "maximal number of separation rounds in the root node (-1: unlimited)",
         &sepadata->maxroundsroot, FALSE, DEFAULT_MAXROUNDSROOT, -1, INT_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddIntParam(scip,
         "separating/flowcover/maxtries",
         "maximal number of rows to separate flow cover cuts for per separation round (-1: unlimited)",
         &sepadata->maxtries, TRUE, DEFAULT_MAXTRIES, -1, INT_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddIntParam(scip,
         "separating/flowcover/maxtriesroot",
         "maximal number of rows to separate flow cover cuts for per separation round in the root (-1: unlimited)",
         &sepadata->maxtriesroot, TRUE, DEFAULT_MAXTRIESROOT, -1, INT_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddIntParam(scip,
         "separating/flowcover/maxfails",
         "maximal number of consecutive fails to generate a cut per separation round (-1: unlimited)",
         &sepadata->maxfails, TRUE, DEFAULT_MAXFAILS, -1, INT_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddIntParam(scip,
         "separating/flowcover/maxfailsroot",
         "maximal number of consecutive fails to generate a cut per separation round in the root (-1: unlimited)",
         &sepadata->maxfailsroot, TRUE, DEFAULT_MAXFAILSROOT, -1, INT_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddIntParam(scip,
         "separating/flowcover/maxsepacuts",
         "maximal number of flow cover cuts separated per separation round",
         &sepadata->maxsepacuts, FALSE, DEFAULT_MAXSEPACUTS, 0, INT_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddIntParam(scip,
         "separating/flowcover/maxsepacutsroot",
         "maximal number of flow cover cuts separated per separation round in the root",
         &sepadata->maxsepacutsroot, FALSE, DEFAULT_MAXSEPACUTSROOT, 0, INT_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddRealParam(scip,
         "separating/flowcover/maxslack",
         "maximal slack of rows to separate flow cover cuts for",
         &sepadata->maxslack, TRUE, DEFAULT_MAXSLACK, 0.0, SCIP_REAL_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddRealParam(scip,
         "separating/flowcover/maxslackroot",
         "maximal slack of rows to separate flow cover cuts for in the root",
         &sepadata->maxslackroot, TRUE, DEFAULT_MAXSLACKROOT, 0.0, SCIP_REAL_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddRealParam(scip,
         "separating/flowcover/slackscore",
         "weight of slack in the scoring of the rows",
         &sepadata->slackscore, TRUE, DEFAULT_SLACKSCORE, 0.0, SCIP_REAL_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddRealParam(scip,
         "separating/flowcover/maxrowdensity",
         "maximal density of row to separate flow cover cuts for",
         &sepadata->maxrowdensity, TRUE, DEFAULT_MAXROWDENSITY, 0.0, 1.0, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip,
         "separating/flowcover/dynamiccuts",
         "should generated cuts be removed from the LP if they are no longer tight?",
         &sepadata->dynamiccuts, FALSE, DEFAULT_DYNAMICCUTS, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip,
         "separating/flowcover/multbyminusone",
         "should flow cover cuts be separated for 0-1 single node flow set with reversed arcs in addition?",
         &sepadata->multbyminusone, TRUE, DEFAULT_MULTBYMINUSONE, NULL, NULL) );
   SCIP_CALL( SCIPaddIntParam(scip,
         "separating/flowcover/maxtestdelta",
         "cut generation heuristic: maximal number of different deltas to try",
         &sepadata->maxtestdelta, TRUE, DEFAULT_MAXTESTDELTA, 0, INT_MAX, NULL, NULL) );
   return SCIP_OKAY;
}
