/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2015 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   presol_tworowbnd.c
 * @brief  do bound tightening by using two rows
 * @author Dieter Weninger
 *
 * Perform bound tightening on two inequalities with some common variables.
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <stdio.h>
#include <assert.h>
#include <string.h>

#include "scip/scipdefplugins.h"
#include "scip/pub_matrix.h"

#include "presol_tworowbnd.h"

#define PRESOL_NAME            "tworowbnd"
#define PRESOL_DESC            "do bound tigthening by using two rows"
#define PRESOL_PRIORITY           500000     /**< priority of the presolver (>= 0: before, < 0: after constraint handlers) */
#define PRESOL_MAXROUNDS              -1     /**< maximal number of presolving rounds the presolver participates in (-1: no limit) */
#define PRESOL_TIMING           SCIP_PRESOLTIMING_MEDIUM /* timing of the presolver (fast, medium, or exhaustive) */

#define SUPPORT_THRESHOLD            0.5     /**< threshold for two constraints overlap */
#define FASTMODE_THRESHOLD          1000     /**< max number of baserows for switching to fast mode */


/** type of bound change */
enum Bndchgtype
{
   NOCHANGE   = 0,
   LOWERBOUND = 1,
   UPPERBOUND = 2,
   BOTHBOUNDS = 3
};
typedef enum Bndchgtype BNDCHGTYPE;

/*
 * Local methods
 */

#if 0
/** write min and max LP to file */
static
void writeLPs(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIPMILPMATRIX*       matrix,             /**< constraint matrix object */
   int                   baserow,            /**< base row index */
   int                   otherrow,           /**< other row index */
   int                   numoverlap,         /**< overlap-size */
   int*                  overlapidx,         /**< overlap column indexes */
   int*                  othernonoverlapidx, /**< other row non overlap indexes */
   int*                  basenonoverlapidx,  /**< base row non overlap indexes */
   SCIP_Real*            coebaseoverlap,     /**< base row overlap coefficients */
   SCIP_Real*            coeotheroverlap,    /**< other row overlap coefficients */
   SCIP_Real*            coebasenonoverlap,  /**< base row non overlap coefficients */
   SCIP_Real*            coeothernonoverlap, /**< other row non overlap coefficients */
   SCIP_Real*            lowerbds,           /**< lower bounds */
   SCIP_Real*            upperbds            /**< upper bounds */
   )
{
   FILE* filemax;
   FILE* filemin;
   SCIP_Real lhs;
   int i;
   int nothernonolap;

   lhs = SCIPmatrixGetRowLhs(matrix, otherrow);

   filemax = fopen("max.lp", "wt");
   filemin = fopen("min.lp", "wt");
   if( filemax != NULL && filemin != NULL )
   {
      fprintf(filemax,"max\n\t");
      fprintf(filemin,"min\n\t");

      for(i = 0; i < numoverlap; i++)
      {
         if(coebaseoverlap[i] > 0.0)
         {
            fprintf(filemax,"+%f %s ",coebaseoverlap[i],SCIPvarGetName(SCIPmatrixGetVar(matrix,overlapidx[i])));
            fprintf(filemin,"+%f %s ",coebaseoverlap[i],SCIPvarGetName(SCIPmatrixGetVar(matrix,overlapidx[i])));
         }
         else
         {
            fprintf(filemax,"%f %s ",coebaseoverlap[i],SCIPvarGetName(SCIPmatrixGetVar(matrix,overlapidx[i])));
            fprintf(filemin,"%f %s ",coebaseoverlap[i],SCIPvarGetName(SCIPmatrixGetVar(matrix,overlapidx[i])));
         }
      }

      fprintf(filemax,"\ns.t.\n\t");
      fprintf(filemin,"\ns.t.\n\t");

      for(i = 0; i < numoverlap; i++)
      {
         if(coeotheroverlap[i] > 0.0)
         {
            fprintf(filemax,"+%f %s ",coeotheroverlap[i],SCIPvarGetName(SCIPmatrixGetVar(matrix,overlapidx[i])));
            fprintf(filemin,"+%f %s ",coeotheroverlap[i],SCIPvarGetName(SCIPmatrixGetVar(matrix,overlapidx[i])));
         }
         else
         {
            fprintf(filemax,"%f %s ",coeotheroverlap[i],SCIPvarGetName(SCIPmatrixGetVar(matrix,overlapidx[i])));
            fprintf(filemin,"%f %s ",coeotheroverlap[i],SCIPvarGetName(SCIPmatrixGetVar(matrix,overlapidx[i])));
         }
      }

      nothernonolap = SCIPmatrixGetRowNNonzs(matrix, otherrow) - numoverlap;

      for(i = 0; i < nothernonolap; i++)
      {
         if(coeothernonoverlap[i] > 0.0)
         {
            fprintf(filemax,"+%f %s ",coeothernonoverlap[i],SCIPvarGetName(SCIPmatrixGetVar(matrix,othernonoverlapidx[i])));
            fprintf(filemin,"+%f %s ",coeothernonoverlap[i],SCIPvarGetName(SCIPmatrixGetVar(matrix,othernonoverlapidx[i])));
         }
         else
         {
            fprintf(filemax,"%f %s ",coeothernonoverlap[i],SCIPvarGetName(SCIPmatrixGetVar(matrix,othernonoverlapidx[i])));
            fprintf(filemin,"%f %s ",coeothernonoverlap[i],SCIPvarGetName(SCIPmatrixGetVar(matrix,othernonoverlapidx[i])));
         }
      }
      fprintf(filemax," >= %f\n",lhs);
      fprintf(filemin," >= %f\n",lhs);

      fprintf(filemax,"bounds\n");
      fprintf(filemin,"bounds\n");

      for(i = 0; i < numoverlap; i++)
      {
         if(!SCIPisInfinity(scip,-lowerbds[overlapidx[i]]) && !SCIPisInfinity(scip,upperbds[overlapidx[i]]))
         {
            fprintf(filemax,"\t%f <= %s <= %f\n",lowerbds[overlapidx[i]],
               SCIPvarGetName(SCIPmatrixGetVar(matrix,overlapidx[i])),upperbds[overlapidx[i]]);
            fprintf(filemin,"\t%f <= %s <= %f\n",lowerbds[overlapidx[i]],
               SCIPvarGetName(SCIPmatrixGetVar(matrix,overlapidx[i])),upperbds[overlapidx[i]]);
         }
         else if(!SCIPisInfinity(scip,-lowerbds[overlapidx[i]]) )
         {
            fprintf(filemax,"\t%f <= %s\n",lowerbds[overlapidx[i]],SCIPvarGetName(SCIPmatrixGetVar(matrix,overlapidx[i])));
            fprintf(filemin,"\t%f <= %s\n",lowerbds[overlapidx[i]],SCIPvarGetName(SCIPmatrixGetVar(matrix,overlapidx[i])));
         }
         else if(!SCIPisInfinity(scip,upperbds[overlapidx[i]]) )
         {
            fprintf(filemax,"\t%s <= %f\n",SCIPvarGetName(SCIPmatrixGetVar(matrix,overlapidx[i])),upperbds[overlapidx[i]]);
            fprintf(filemin,"\t%s <= %f\n",SCIPvarGetName(SCIPmatrixGetVar(matrix,overlapidx[i])),upperbds[overlapidx[i]]);
         }
      }

      for(i = 0; i < nothernonolap; i++)
      {
         if(!SCIPisInfinity(scip,-lowerbds[othernonoverlapidx[i]]) && !SCIPisInfinity(scip,upperbds[othernonoverlapidx[i]]))
         {
            fprintf(filemax,"\t%f <= %s <= %f\n",lowerbds[othernonoverlapidx[i]],
               SCIPvarGetName(SCIPmatrixGetVar(matrix,othernonoverlapidx[i])),upperbds[othernonoverlapidx[i]]);
            fprintf(filemin,"\t%f <= %s <= %f\n",lowerbds[othernonoverlapidx[i]],
               SCIPvarGetName(SCIPmatrixGetVar(matrix,othernonoverlapidx[i])),upperbds[othernonoverlapidx[i]]);
         }
         else if(!SCIPisInfinity(scip,-lowerbds[othernonoverlapidx[i]]) )
         {
            fprintf(filemax,"\t%f <= %s\n",lowerbds[othernonoverlapidx[i]],
               SCIPvarGetName(SCIPmatrixGetVar(matrix,othernonoverlapidx[i])));
            fprintf(filemin,"\t%f <= %s\n",lowerbds[othernonoverlapidx[i]],
               SCIPvarGetName(SCIPmatrixGetVar(matrix,othernonoverlapidx[i])));
         }
         else if(!SCIPisInfinity(scip,upperbds[othernonoverlapidx[i]]) )
         {
            fprintf(filemax,"\t%s <= %f\n",SCIPvarGetName(SCIPmatrixGetVar(matrix,othernonoverlapidx[i])),upperbds[othernonoverlapidx[i]]);
            fprintf(filemin,"\t%s <= %f\n",SCIPvarGetName(SCIPmatrixGetVar(matrix,othernonoverlapidx[i])),upperbds[othernonoverlapidx[i]]);
         }
      }


      fprintf(filemax,"end\n");
      fprintf(filemin,"end\n");

      fclose(filemax);
      fclose(filemin);
   }
   else
      assert(0);
}
#endif


/** solve two LPs with one row each
 *
 * a1x + a3y      >= b1  (other row)
 * a2x      + a4z >= b2  (base row)
 *
 * minact = min{a2x : a1x + a3y >= b1}
 * maxact = max{a2x : a1x + a3y >= b1}
 */
static
void getactivities(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIPMILPMATRIX*       matrix,             /**< constraint matrix object */
   int                   baserow,            /**< base row index */
   int                   otherrow,           /**< other row index */
   int                   numoverlap,         /**< overlap-size */
   int*                  overlapidx,         /**< overlap column indexes */
   int*                  othernonoverlapidx, /**< other row non overlap indexes */
   int*                  basenonoverlapidx,  /**< base row non overlap indexes */
   SCIP_Real*            coebaseoverlap,     /**< base row overlap coefficients */
   SCIP_Real*            coeotheroverlap,    /**< other row overlap coefficients */
   SCIP_Real*            coebasenonoverlap,  /**< base row non overlap coefficients */
   SCIP_Real*            coeothernonoverlap, /**< other row non overlap coefficients */
   SCIP_Real*            lowerbds,           /**< lower bounds */
   SCIP_Real*            upperbds,           /**< upper bounds */
   SCIP_Real*            tmplowerbds,        /**< tmp lower bounds */
   SCIP_Real*            tmpupperbds,        /**< tmp upper bounds */
   SCIP_Real*            minratios,          /**< min LP ratios */
   SCIP_Real*            maxratios,          /**< max LP ratios */
   int*                  minsortedidx,       /**< min LP sorted indexes */
   int*                  maxsortedidx,       /**< max LP sorted indexes */
   SCIP_Real*            minact,             /**< calculated overlap minimal activity w.r.t. to the other row */
   SCIP_Real*            maxact              /**< calculated overlap maximal activity w.r.t. to the other row */
   )
{
   SCIP_Real val;
   int nothernonolap;
   SCIP_Real lhs;
   SCIP_Real minlhs;
   SCIP_Real maxlhs;
   SCIP_Bool consred;
   SCIP_Real objoffset;
   int nminratios;
   int nmaxratios;
   int i;
   SCIP_Bool infpresent;

   *minact = 0;
   *maxact = 0;
   nminratios = 0;
   nmaxratios = 0;

#if 0
   writeLPs(scip,matrix,baserow,otherrow,numoverlap,overlapidx,othernonoverlapidx,basenonoverlapidx,
      coebaseoverlap,coeotheroverlap,coebasenonoverlap,coeothernonoverlap,lowerbds,upperbds);
#endif

   lhs = SCIPmatrixGetRowLhs(matrix, otherrow);
   assert(lhs > -SCIPinfinity(scip));

   nothernonolap = SCIPmatrixGetRowNNonzs(matrix, otherrow) - numoverlap;
   val = 0;
   infpresent = FALSE;
   consred = FALSE;

   for(i = 0; i < nothernonolap; i++)
   {
      if(coeothernonoverlap[i] < 0.0)
      {
         if( SCIPisInfinity(scip, -lowerbds[othernonoverlapidx[i]]) )
         {
            consred = TRUE;
            break;
         }
         else
         {
            val += coeothernonoverlap[i] * lowerbds[othernonoverlapidx[i]];
         }
      }
      else if(coeothernonoverlap[i] > 0.0)
      {
         if( SCIPisInfinity(scip, upperbds[othernonoverlapidx[i]]) )
         {
            consred = TRUE;
            break;
         }
         else
         {
            val += coeothernonoverlap[i] * upperbds[othernonoverlapidx[i]];
         }
      }
   }

   if(!consred)
   {
      lhs -= val;
      objoffset = 0;

      for(i = 0; i < numoverlap; i++)
      {
         tmplowerbds[i] = lowerbds[overlapidx[i]];
         tmpupperbds[i] = upperbds[overlapidx[i]];

         if(coeotheroverlap[i] > 0.0)
         {
            SCIP_Real tmp;
            tmp = tmplowerbds[i];
            tmplowerbds[i] = -tmpupperbds[i];
            tmpupperbds[i] = -tmp;

            coeotheroverlap[i] = -coeotheroverlap[i];
            coebaseoverlap[i] = -coebaseoverlap[i];
         }

         if(coebaseoverlap[i] < 0.0)
         {
            minratios[nminratios] = coebaseoverlap[i] / coeotheroverlap[i];
            minsortedidx[nminratios] = i;
            nminratios++;
         }
         else
         {
            maxratios[nmaxratios] = coebaseoverlap[i] / coeotheroverlap[i];
            maxsortedidx[nmaxratios] = i;
            nmaxratios++;
         }

         if( SCIPisInfinity(scip, -tmplowerbds[i]) || SCIPisInfinity(scip, tmplowerbds[i]) )
         {
            infpresent = TRUE;
            *minact = -SCIPinfinity(scip);
            *maxact = SCIPinfinity(scip);
            break;
         }
         else
         {
            lhs -= coeotheroverlap[i] * tmplowerbds[i];
            objoffset += coebaseoverlap[i] * tmplowerbds[i];
            tmpupperbds[i] -= tmplowerbds[i];
         }
      }

      if( !infpresent )
      {
         /* min case */
         SCIPsortRealInt(minratios, minsortedidx, nminratios);
         minlhs = lhs;
         for( i = nminratios-1; 0 <= i; i-- )
         {
            SCIP_Real tmpval;

            if( SCIPisInfinity(scip, tmpupperbds[minsortedidx[i]]) )
            {
               infpresent = TRUE;
               *minact = -SCIPinfinity(scip);
               break;
            }
            else
            {
               tmpval = coeotheroverlap[minsortedidx[i]] * tmpupperbds[minsortedidx[i]];
               if( SCIPisGE(scip, tmpval, minlhs) )
               {
                  *minact += coebaseoverlap[minsortedidx[i]] * tmpupperbds[minsortedidx[i]];
                  minlhs -= tmpval;
               }
               else
               {
                  tmpval = minlhs / coeotheroverlap[minsortedidx[i]];
                  if(tmpval > 0.0)
                  {
                     assert(tmpval <= tmpupperbds[minsortedidx[i]]);
                     *minact += coebaseoverlap[minsortedidx[i]] * tmpval;
                  }
                  break;
               }
            }
         }
         if( !infpresent )
            *minact += objoffset;

         /* max case */
         SCIPsortRealInt(maxratios, maxsortedidx, nmaxratios);
         maxlhs = lhs;
         for( i = 0; i < nmaxratios; i++ )
         {
            SCIP_Real tmpval;

            if( SCIPisInfinity(scip, tmpupperbds[maxsortedidx[i]]) )
            {
               infpresent = TRUE;
               *maxact = SCIPinfinity(scip);
               break;
            }
            else
            {
               tmpval = coeotheroverlap[maxsortedidx[i]] * tmpupperbds[maxsortedidx[i]];
               if( SCIPisGE(scip, tmpval, maxlhs) )
               {
                  *maxact += coebaseoverlap[maxsortedidx[i]] * tmpupperbds[maxsortedidx[i]];
                  maxlhs -= tmpval;
               }
               else
               {
                  tmpval = maxlhs / coeotheroverlap[maxsortedidx[i]];
                  if(tmpval > 0.0)
                  {
                     assert(tmpval <= tmpupperbds[maxsortedidx[i]]);
                     *maxact += coebaseoverlap[maxsortedidx[i]] * tmpval;
                  }
                  break;
               }
            }
         }
         if( !infpresent )
            *maxact += objoffset;
      }
   }
   else
   {
      /* min case */
      for(i = 0; i < numoverlap; i++)
      {
         if(coebaseoverlap[i] > 0.0)
         {
            if( !SCIPisInfinity(scip, -lowerbds[overlapidx[i]]) )
            {
               *minact +=  coebaseoverlap[i] * lowerbds[overlapidx[i]];
            }
            else
            {
               *minact = -SCIPinfinity(scip);
               break;
            }
         }
         else if(coebaseoverlap[i] < 0.0)
         {
            if( !SCIPisInfinity(scip, upperbds[overlapidx[i]]) )
            {
               *minact +=  coebaseoverlap[i] * upperbds[overlapidx[i]];
            }
            else
            {
               *minact = -SCIPinfinity(scip);
               break;
            }
         }
      }

      /* max case */
      for(i = 0; i < numoverlap; i++)
      {
         if(coebaseoverlap[i] > 0.0)
         {
            if( !SCIPisInfinity(scip, upperbds[overlapidx[i]]) )
            {
               *maxact +=  coebaseoverlap[i] * upperbds[overlapidx[i]];
            }
            else
            {
               *maxact = SCIPinfinity(scip);
               break;
            }
         }
         else if(coebaseoverlap[i] < 0.0)
         {
            if( !SCIPisInfinity(scip, -lowerbds[overlapidx[i]]) )
            {
               *maxact +=  coebaseoverlap[i] * lowerbds[overlapidx[i]];
            }
            else
            {
               *maxact = SCIPinfinity(scip);
               break;
            }
         }
      }
   }
}

/** calculate min activity */
static
SCIP_Real getinfimum(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIPMILPMATRIX*       matrix,             /**< constraint matrix object */
   int                   len,                /**< length */
   int*                  varidxs,            /**< variables indexes */
   SCIP_Real*            coeffs,             /**< coefficients */
   SCIP_Real*            lowerbds,           /**< lower bounds */
   SCIP_Real*            upperbds            /**< upper bounds */
   )
{
   SCIP_Real infimum;
   int i;

   infimum = 0;

   for(i = 0; i < len; i++)
   {
      if(coeffs[i] > 0.0)
      {
         if( SCIPisInfinity(scip,-lowerbds[varidxs[i]]) )
         {
            infimum = -SCIPinfinity(scip);
            break;
         }
         else
         {
            infimum += coeffs[i] * lowerbds[varidxs[i]];
         }
      }
      else
      {
         if( SCIPisInfinity(scip,upperbds[varidxs[i]]) )
         {
            infimum = -SCIPinfinity(scip);
            break;
         }
         else
         {
            infimum += coeffs[i] * upperbds[varidxs[i]];
         }
      }
   }

   return infimum;
}

/**< calculate max activity */
static
SCIP_Real getsupremum(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIPMILPMATRIX*       matrix,             /**< constraint matrix object */
   int                   len,                /**< length */
   int*                  varidxs,            /**< variable indexes */
   SCIP_Real*            coeffs,             /**< coefficients */
   SCIP_Real*            lowerbds,           /**< lower bounds */
   SCIP_Real*            upperbds,           /**< upper bounds */
   int*                  infcnt              /**< infinity counter */
   )
{
   SCIP_Real supremum;
   int i;

   *infcnt = 0;
   supremum = 0;

   for(i = 0; i < len; i++)
   {
      if(coeffs[i] < 0.0)
      {
         if( SCIPisInfinity(scip,-lowerbds[varidxs[i]]) )
            (*infcnt)++;
         else
            supremum += coeffs[i] * lowerbds[varidxs[i]];
      }
      else
      {
         if( SCIPisInfinity(scip,upperbds[varidxs[i]]) )
            (*infcnt)++;
         else
            supremum += coeffs[i] * upperbds[varidxs[i]];
      }
   }

   if(*infcnt > 0)
      supremum = SCIPinfinity(scip);

   return supremum;
}

/**< get max activity without one column */
static
SCIP_Real getsupremumidx(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIPMILPMATRIX*       matrix,             /**< constraint matrix object */
   int                   len,                /**< length */
   int*                  varidxs,            /**< variable indexes */
   SCIP_Real*            coeffs,             /**< coefficients */
   SCIP_Real*            lowerbds,           /**< upper bounds */
   SCIP_Real*            upperbds,           /**< lower bounds */
   int                   idx                 /**< omitting index */
   )
{
   SCIP_Real supremum;
   int i;

   supremum = 0;

   for(i = 0; i < len; i++)
   {
      if(i == idx)
         continue;

      if(coeffs[i] < 0.0)
      {
         assert(!SCIPisInfinity(scip,-lowerbds[varidxs[i]]));
         supremum += coeffs[i] * lowerbds[varidxs[i]];
      }
      else
      {
         assert(!SCIPisInfinity(scip,upperbds[varidxs[i]]));
         supremum += coeffs[i] * upperbds[varidxs[i]];
      }
   }

   return supremum;
}

/** apply bound tightening on two overlapping constraints */
static
void applytightening(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIPMILPMATRIX*       matrix,             /**< constraint matrix object */
   int                   baserow,            /**< base row index */
   int                   otherrow,           /**< other row index */
   int                   numoverlap,         /**< overlap-size */
   int*                  overlapidx,         /**< overlap column indexes */
   int*                  othernonoverlapidx, /**< other row non overlap indexes */
   int*                  basenonoverlapidx,  /**< base row non overlap indexes */
   SCIP_Real*            coebaseoverlap,     /**< base row overlap coefficients */
   SCIP_Real*            coeotheroverlap,    /**< other row overlap coefficients */
   SCIP_Real*            coebasenonoverlap,  /**< base row non overlap coefficients */
   SCIP_Real*            coeothernonoverlap, /**< other row non overlap coefficients */
   SCIP_Real*            lowerbds,           /**< lower bounds */
   SCIP_Real*            upperbds,           /**< upper bounds */
   SCIP_Real*            tmplowerbds,        /**< tmp lower bounds */
   SCIP_Real*            tmpupperbds,        /**< tmp upper bounds */
   SCIP_Real*            minratios,          /**< min LP ratios */
   SCIP_Real*            maxratios,          /**< max LP ratios */
   int*                  minsortedidx,       /**< min LP sorted indexes */
   int*                  maxsortedidx,       /**< max LP sorted indexes */
   int*                  ntightenbnds,       /**< number of tightened bounds */
   BNDCHGTYPE*           tighten,            /**< tightened bounds */
   int*                  ndeletecons,        /**< number of redundant constraints */
   SCIP_Bool*            deletecons          /**< redundant constraints */
   )
{
   SCIP_Real maxact;
   SCIP_Real minact;
   SCIP_Real infimum;
   SCIP_Real supremum;
   int len;
   SCIP_Real lhs;
   int i;

   getactivities(scip, matrix, baserow, otherrow, numoverlap, overlapidx,
      othernonoverlapidx, basenonoverlapidx,
      coebaseoverlap, coeotheroverlap, coebasenonoverlap, coeothernonoverlap,
      lowerbds, upperbds, tmplowerbds, tmpupperbds, minratios, maxratios,
      minsortedidx, maxsortedidx, &minact, &maxact);

   len = SCIPmatrixGetRowNNonzs(matrix, baserow) - numoverlap;
   lhs = SCIPmatrixGetRowLhs(matrix, baserow);

   if( !SCIPisInfinity(scip, -minact) )
   {
      /* detect redundant constraints */
      infimum = getinfimum(scip, matrix, len, basenonoverlapidx, coebasenonoverlap, lowerbds, upperbds);
      if( !SCIPisInfinity(scip, -infimum) )
      {
         if( SCIPisGE(scip, minact+infimum, lhs) )
         {
            if( !deletecons[baserow] )
            {
               (*ndeletecons)++;
               deletecons[baserow] = TRUE;
            }
         }
      }
   }

   if( !SCIPisInfinity(scip, maxact) )
   {
      int infcnt;
      SCIP_Real bnd;
      SCIP_Real tmpsup;

      /* bound tightening */
      supremum = getsupremum(scip, matrix, len, basenonoverlapidx, coebasenonoverlap, lowerbds, upperbds, &infcnt);
      if( !SCIPisInfinity(scip, supremum) )
      {
         for(i = 0; i < len; i++)
         {
            if(coebasenonoverlap[i] < 0.0)
            {
               /* get ub */
               tmpsup = supremum - (coebasenonoverlap[i] * lowerbds[basenonoverlapidx[i]]);
               bnd = (lhs - (tmpsup + maxact)) / coebasenonoverlap[i];
               if(bnd < upperbds[basenonoverlapidx[i]])
               {
                  upperbds[basenonoverlapidx[i]] = bnd;
                  if(tighten[basenonoverlapidx[i]] != UPPERBOUND && tighten[basenonoverlapidx[i]] != BOTHBOUNDS)
                  {
                     (*ntightenbnds)++;
                     if(tighten[basenonoverlapidx[i]] == LOWERBOUND)
                        tighten[basenonoverlapidx[i]] = BOTHBOUNDS;
                     else
                        tighten[basenonoverlapidx[i]] = UPPERBOUND;
                  }
               }
            }
            else
            {
               /* get lb */
               tmpsup = supremum - (coebasenonoverlap[i] * upperbds[basenonoverlapidx[i]]);
               bnd = (lhs - (tmpsup + maxact)) / coebasenonoverlap[i];
               if(bnd > lowerbds[basenonoverlapidx[i]])
               {
                  lowerbds[basenonoverlapidx[i]] = bnd;
                  if(tighten[basenonoverlapidx[i]] != LOWERBOUND && tighten[basenonoverlapidx[i]] != BOTHBOUNDS)
                  {
                     (*ntightenbnds)++;
                     if(tighten[basenonoverlapidx[i]] == UPPERBOUND)
                        tighten[basenonoverlapidx[i]] = BOTHBOUNDS;
                     else
                        tighten[basenonoverlapidx[i]] = LOWERBOUND;
                  }
               }
            }
         }
      }
      else
      {
         for(i = 0; i < len; i++)
         {
            if(coebasenonoverlap[i] < 0.0)
            {
               if(infcnt == 1 && SCIPisInfinity(scip, -lowerbds[basenonoverlapidx[i]]) )
               {
                  /* get ub */
                  tmpsup = getsupremumidx(scip, matrix, len, basenonoverlapidx, coebasenonoverlap, lowerbds, upperbds, i);
                  assert(!SCIPisInfinity(scip, tmpsup));
                  bnd = (lhs - (tmpsup + maxact)) / coebasenonoverlap[i];
                  if(bnd < upperbds[basenonoverlapidx[i]])
                  {
                     upperbds[basenonoverlapidx[i]] = bnd;
                     if(tighten[basenonoverlapidx[i]] != UPPERBOUND && tighten[basenonoverlapidx[i]] != BOTHBOUNDS)
                     {
                        (*ntightenbnds)++;
                        if(tighten[basenonoverlapidx[i]] == LOWERBOUND)
                           tighten[basenonoverlapidx[i]] = BOTHBOUNDS;
                        else
                           tighten[basenonoverlapidx[i]] = UPPERBOUND;
                     }
                  }
               }
            }
            else
            {
               if(infcnt == 1 && SCIPisInfinity(scip, upperbds[basenonoverlapidx[i]]) )
               {
                  /* get lb */
                  tmpsup = getsupremumidx(scip, matrix, len, basenonoverlapidx, coebasenonoverlap, lowerbds, upperbds, i);
                  assert(!SCIPisInfinity(scip, tmpsup));
                  bnd = (lhs - (tmpsup + maxact)) / coebasenonoverlap[i];
                  if(bnd > lowerbds[basenonoverlapidx[i]])
                  {
                     lowerbds[basenonoverlapidx[i]] = bnd;
                     if(tighten[basenonoverlapidx[i]] != LOWERBOUND && tighten[basenonoverlapidx[i]] != BOTHBOUNDS)
                     {
                        (*ntightenbnds)++;
                        if(tighten[basenonoverlapidx[i]] == UPPERBOUND)
                           tighten[basenonoverlapidx[i]] = BOTHBOUNDS;
                        else
                           tighten[basenonoverlapidx[i]] = LOWERBOUND;
                     }
                  }
               }
            }
         }
      }
   }
}

/** extract coefficients from matrix */
static
void getcoefficients(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIPMILPMATRIX*       matrix,             /**< constraint matrix object */
   int                   baserow,            /**< base row index */
   int                   otherrow,           /**< other row index */
   int                   numoverlap,         /**< overlap-size */
   int*                  olapidxbaseorder,   /**< overlap column indexes in baserow order */
   int*                  olapidxotherorder,  /**< overlap column indexes in otherrow order */
   int*                  othernonoverlapidx, /**< other row non overlap indexes */
   int*                  basenonoverlapidx,  /**< base row non overlap indexes */
   SCIP_Real*            coebaseoverlap,     /**< base row overlap coefficients */
   SCIP_Real*            coeotheroverlap,    /**< other row overlap coefficients */
   SCIP_Real*            coebasenonoverlap,  /**< base row non overlap coefficients */
   SCIP_Real*            coeothernonoverlap  /**< other row non overlap coefficients */
   )
{
   SCIP_Real* valpnt;
   int* rowpnt;
   int* rowend;
   int baserowcnt;
   int otherrowcnt;
   int olapcnt;
   int nonolapcnt;

   baserowcnt = SCIPmatrixGetRowNNonzs(matrix, baserow);
   otherrowcnt = SCIPmatrixGetRowNNonzs(matrix, otherrow);
   assert( baserowcnt != 0 && otherrowcnt != 0 );

   /* set end marker */
   if(numoverlap < SCIPmatrixGetNColumns(matrix))
   {
      olapidxbaseorder[numoverlap] = -1;
      olapidxotherorder[numoverlap] = -1;
   }

   olapcnt = 0;
   nonolapcnt = 0;

   rowpnt = SCIPmatrixGetRowIdxPtr(matrix, baserow);
   rowend = rowpnt + baserowcnt;
   valpnt = SCIPmatrixGetRowValPtr(matrix, baserow);
   for( ; rowpnt < rowend; rowpnt++, valpnt++ )
   {
      if( olapidxbaseorder[olapcnt] == *rowpnt )
      {
         coebaseoverlap[olapcnt] = *valpnt;
         olapcnt++;
      }
      else
      {
         basenonoverlapidx[nonolapcnt] = *rowpnt;
         coebasenonoverlap[nonolapcnt] = *valpnt;
         nonolapcnt++;
      }
   }

   assert(olapcnt+nonolapcnt == baserowcnt);
   assert(olapcnt == numoverlap);
   assert(nonolapcnt > 0);

   olapcnt = 0;
   nonolapcnt = 0;

   rowpnt = SCIPmatrixGetRowIdxPtr(matrix, otherrow);
   rowend = rowpnt + otherrowcnt;
   valpnt = SCIPmatrixGetRowValPtr(matrix, otherrow);
   for( ; rowpnt < rowend; rowpnt++, valpnt++ )
   {
      if( olapidxotherorder[olapcnt] == *rowpnt )
      {
         coeotheroverlap[olapcnt] = *valpnt;
         olapcnt++;
      }
      else
      {
         othernonoverlapidx[nonolapcnt] = *rowpnt;
         coeothernonoverlap[nonolapcnt] = *valpnt;
         nonolapcnt++;
      }
   }

   assert(olapcnt+nonolapcnt == otherrowcnt);
   assert(olapcnt == numoverlap);
}

/** calculate overlap-size */
static
void getnumoverlap(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIPMILPMATRIX*       matrix,             /**< constraint matrix object */
   int                   baserow,            /**< base row index */
   int                   otherrow,           /**< other row index */
   int*                  countings,          /**< overlap counting helper array */
   int*                  clearinfo,          /**< reset helper array */
   int*                  numoverlap,         /**< overlap-size */
   int*                  olapidxotherorder   /**< overlap column indexes in otherrow order */
   )
{
   int* rowpnt;
   int* rowend;
   int noverlap;
   int baserowcnt;
   int otherrowcnt;
   int nclear;
   int i;

   noverlap = 0;
   nclear = 0;

   baserowcnt = SCIPmatrixGetRowNNonzs(matrix, baserow);
   otherrowcnt = SCIPmatrixGetRowNNonzs(matrix, otherrow);
   if( baserowcnt == 0 || otherrowcnt == 0 )
   {
      *numoverlap = noverlap;
      return;
   }

   /* set flags corresponding to baserow non-zeros */
   rowpnt = SCIPmatrixGetRowIdxPtr(matrix, baserow);
   rowend = rowpnt + baserowcnt;
   for( ; rowpnt < rowend; rowpnt++ )
   {
      countings[*rowpnt] = 1;
      clearinfo[nclear] = *rowpnt;
      nclear++;
   }

   rowpnt = SCIPmatrixGetRowIdxPtr(matrix, otherrow);
   rowend = rowpnt + otherrowcnt;
   for( ; rowpnt < rowend; rowpnt++ )
   {
      if(countings[*rowpnt] == 1)
      {
         /* collect overlapping indexes in otherrow order */
         olapidxotherorder[noverlap] = *rowpnt;
         noverlap++;
      }
   }

   for( i = 0; i < nclear; i++ )
      countings[clearinfo[i]] = 0;

   *numoverlap = noverlap;
}

static
void getoverlapbaseordered(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIPMILPMATRIX*       matrix,             /**< constraint matrix object */
   int                   baserow,            /**< base row index */
   int                   otherrow,           /**< other row index */
   int*                  countings,          /**< overlap counting helper array */
   int*                  clearinfo,          /**< reset helper array */
   int                   numoverlap,         /**< just calculated overlap-size */
   int*                  olapidxbaseorder    /**< overlap column indexes in baserow order */
   )
{
   int* rowpnt;
   int* rowend;
   int noverlap;
   int baserowcnt;
   int otherrowcnt;
   int nclear;
   int i;

   noverlap = 0;
   nclear = 0;

   baserowcnt = SCIPmatrixGetRowNNonzs(matrix, baserow);
   otherrowcnt = SCIPmatrixGetRowNNonzs(matrix, otherrow);

   /* set flags corresponding to otherrow non-zeros */
   rowpnt = SCIPmatrixGetRowIdxPtr(matrix, otherrow);
   rowend = rowpnt + otherrowcnt;
   for( ; rowpnt < rowend; rowpnt++ )
   {
      countings[*rowpnt] = 1;
      clearinfo[nclear] = *rowpnt;
      nclear++;
   }

   rowpnt = SCIPmatrixGetRowIdxPtr(matrix, baserow);
   rowend = rowpnt + baserowcnt;
   for( ; rowpnt < rowend; rowpnt++ )
   {
      if(countings[*rowpnt] == 1)
      {
         /* collect overlapping indexes in baserow order */
         olapidxbaseorder[noverlap] = *rowpnt;
         noverlap++;
      }
   }

   for( i = 0; i < nclear; i++ )
      countings[clearinfo[i]] = 0;

   assert(noverlap == numoverlap);
}


/** perform bound tightening on two rows with a specific support intersection */
static
SCIP_RETCODE calctworowbnds(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIPMILPMATRIX*       matrix,             /**< constraint matrix object */
   int                   nbaserows,          /**< number of base rows */
   int*                  baserows,           /**< base rows indexes */
   SCIP_Real*            lowerbds,           /**< lower bounds */
   SCIP_Real*            upperbds,           /**< upper bounds */
   int*                  ntightenbnds,       /**< number of tightened bounds */
   BNDCHGTYPE*           tighten,            /**< bound tighten information */
   int*                  ndeletecons,        /**< number of redundant constraints */
   SCIP_Bool*            deletecons          /**< redundant constraints */
   )
{
   int* rowpnt;
   int* rowend;
   int rowcnt;
   int br;
   int col;
   int* colpnt;
   int* colend;
   int colcnt;
   int numoverlap;
   int* olapidxbaseorder;
   int* olapidxotherorder;
   SCIP_Real threshold;
   int rowcnt2;
   int nrows;
   int ncols;
   int* othernonoverlapidx;
   int* basenonoverlapidx;
   SCIP_Real* coebaseoverlap;
   SCIP_Real* coeotheroverlap;
   SCIP_Real* coebasenonoverlap;
   SCIP_Real* coeothernonoverlap;
   int* countings;
   int* clearinfo;
   SCIP_Real* tmplowerbds;
   SCIP_Real* tmpupperbds;
   SCIP_Real* minratios;
   SCIP_Real* maxratios;
   int* minsortedidx;
   int* maxsortedidx;
   SCIP_Bool usefastmode;
   int* ignorerowidx;
   SCIP_Bool* ignorerow;
   int ignorerowcnt;
   int i;

   nrows = SCIPmatrixGetNRows(matrix);
   ncols = SCIPmatrixGetNColumns(matrix);

   SCIP_CALL( SCIPallocBufferArray(scip, &olapidxbaseorder, ncols) );
   SCIP_CALL( SCIPallocBufferArray(scip, &olapidxotherorder, ncols) );
   SCIP_CALL( SCIPallocBufferArray(scip, &othernonoverlapidx, ncols) );
   SCIP_CALL( SCIPallocBufferArray(scip, &basenonoverlapidx, ncols) );
   SCIP_CALL( SCIPallocBufferArray(scip, &coebaseoverlap, ncols) );
   SCIP_CALL( SCIPallocBufferArray(scip, &coeotheroverlap, ncols) );
   SCIP_CALL( SCIPallocBufferArray(scip, &coebasenonoverlap, ncols) );
   SCIP_CALL( SCIPallocBufferArray(scip, &coeothernonoverlap, ncols) );
   SCIP_CALL( SCIPallocBufferArray(scip, &countings, ncols) );
   BMSclearMemoryArray(countings, ncols);
   SCIP_CALL( SCIPallocBufferArray(scip, &clearinfo, ncols) );
   SCIP_CALL( SCIPallocBufferArray(scip, &tmplowerbds, ncols) );
   SCIP_CALL( SCIPallocBufferArray(scip, &tmpupperbds, ncols) );
   SCIP_CALL( SCIPallocBufferArray(scip, &minratios, ncols) );
   SCIP_CALL( SCIPallocBufferArray(scip, &maxratios, ncols) );
   SCIP_CALL( SCIPallocBufferArray(scip, &minsortedidx, ncols) );
   SCIP_CALL( SCIPallocBufferArray(scip, &maxsortedidx, ncols) );
   SCIP_CALL( SCIPallocBufferArray(scip, &ignorerowidx, nrows) );
   SCIP_CALL( SCIPallocBufferArray(scip, &ignorerow, nrows) );
   BMSclearMemoryArray(ignorerow, nrows);

   /* use fast mode if too much base rows are present */
   if(nbaserows > FASTMODE_THRESHOLD)
      usefastmode = TRUE;
   else
      usefastmode = FALSE;

   for( br = 0; br < nbaserows; br++ )
   {
      ignorerowcnt = 0;

      rowpnt = SCIPmatrixGetRowIdxPtr(matrix, baserows[br]);
      rowcnt = SCIPmatrixGetRowNNonzs(matrix, baserows[br]);
      if( rowcnt == 0 )
         continue;

      rowend = rowpnt + rowcnt;

      for(; (rowpnt < rowend); rowpnt++ )
      {
         col = *rowpnt;
         colpnt = SCIPmatrixGetColIdxPtr(matrix, col);
         colcnt = SCIPmatrixGetColNNonzs(matrix, col);
         colend = colpnt + colcnt;
         for(; (colpnt < colend); colpnt++ )
         {
            if( *colpnt == baserows[br] || ignorerow[*colpnt] )
               continue;

            /* we consider only >= constraints */
            if( !SCIPmatrixIsRowRhsInfinity(matrix, *colpnt) )
               continue;

            /* determine overlap-size */
            getnumoverlap(scip, matrix, baserows[br], *colpnt,
               countings, clearinfo, &numoverlap, olapidxotherorder);

            if( numoverlap == 0 )
               continue;

            rowcnt2 = SCIPmatrixGetRowNNonzs(matrix, *colpnt);
            threshold = (double)numoverlap/(double)MIN(rowcnt, rowcnt2);

            /* verify if overlap-size is ok */
            if( SUPPORT_THRESHOLD <= threshold && numoverlap < rowcnt )
            {
               getoverlapbaseordered(scip, matrix, baserows[br], *colpnt,
                  countings, clearinfo, numoverlap, olapidxbaseorder);

               getcoefficients(scip, matrix, baserows[br], *colpnt, numoverlap,
                  olapidxbaseorder, olapidxotherorder, othernonoverlapidx, basenonoverlapidx,
                  coebaseoverlap, coeotheroverlap, coebasenonoverlap, coeothernonoverlap);

               applytightening(scip, matrix, baserows[br], *colpnt, numoverlap, olapidxotherorder,
                  othernonoverlapidx, basenonoverlapidx,
                  coebaseoverlap, coeotheroverlap, coebasenonoverlap, coeothernonoverlap,
                  lowerbds, upperbds, tmplowerbds, tmpupperbds, minratios, maxratios,
                  minsortedidx, maxsortedidx, ntightenbnds, tighten, ndeletecons, deletecons);
            }

            ignorerow[*colpnt] = TRUE;
            ignorerowidx[ignorerowcnt] = *colpnt;
            ignorerowcnt++;

            if(usefastmode)
               break;
         }

         if(usefastmode)
            break;
      }

      for(i = 0; i < ignorerowcnt; i++)
         ignorerow[ignorerowidx[i]] = FALSE;
   }

   SCIPfreeBufferArray(scip, &ignorerow);
   SCIPfreeBufferArray(scip, &ignorerowidx);
   SCIPfreeBufferArray(scip, &maxsortedidx);
   SCIPfreeBufferArray(scip, &minsortedidx);
   SCIPfreeBufferArray(scip, &maxratios);
   SCIPfreeBufferArray(scip, &minratios);
   SCIPfreeBufferArray(scip, &tmpupperbds);
   SCIPfreeBufferArray(scip, &tmplowerbds);
   SCIPfreeBufferArray(scip, &clearinfo);
   SCIPfreeBufferArray(scip, &countings);
   SCIPfreeBufferArray(scip, &coeothernonoverlap);
   SCIPfreeBufferArray(scip, &coebasenonoverlap);
   SCIPfreeBufferArray(scip, &coeotheroverlap);
   SCIPfreeBufferArray(scip, &coebaseoverlap);
   SCIPfreeBufferArray(scip, &basenonoverlapidx);
   SCIPfreeBufferArray(scip, &othernonoverlapidx);
   SCIPfreeBufferArray(scip, &olapidxotherorder);
   SCIPfreeBufferArray(scip, &olapidxbaseorder);

   return SCIP_OKAY;
}

/** determine base rows */
static
SCIP_RETCODE getBaseRows(
   SCIP*                 scip,               /**< SCIP main data structure */
   SCIPMILPMATRIX*       matrix,             /**< constraint matrix */
   int*                  nbaserows,          /**< number of present base rows */
   int*                  baserows            /**< indexes of base rows */
   )
{
   int nrows;
   int fill;
   int r;

   assert(scip != NULL);
   assert(matrix != NULL);
   assert(nbaserows != NULL);
   assert(baserows != NULL);

   nrows = SCIPmatrixGetNRows(matrix);

   fill = 0;
   for( r = 0; r < nrows; r++ )
   {
      if( !SCIPmatrixIsRowRhsInfinity(matrix, r) )
         continue;

      baserows[fill] = r;
      fill++;
   }

   *nbaserows = fill;

   return SCIP_OKAY;
}

/** get bounds of variables */
static
void getBounds(
   SCIP*                 scip,               /**< SCIP main data structure */
   SCIPMILPMATRIX*       matrix,             /**< constraint matrix */
   SCIP_Real*            lowerbds,           /**< lower bounds */
   SCIP_Real*            upperbds            /**< upper bounds */
   )
{
   int c;
   int ncols;

   assert(scip != NULL);
   assert(matrix != NULL);
   assert(lowerbds != NULL);
   assert(upperbds != NULL);

   ncols = SCIPmatrixGetNColumns(matrix);

   for( c = 0; c < ncols; c++ )
   {
      SCIP_VAR* var;
      var = SCIPmatrixGetVar(matrix, c);
      lowerbds[c] = SCIPvarGetLbGlobal(var);
      upperbds[c] = SCIPvarGetUbGlobal(var);
   }
}

/*
 * Callback methods of presolver
 */

/** copy method for constraint handler plugins (called when SCIP copies plugins) */
static
SCIP_DECL_PRESOLCOPY(presolCopyTworowbnd)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(presol != NULL);
   assert(strcmp(SCIPpresolGetName(presol), PRESOL_NAME) == 0);

   /* call inclusion method of presolver */
   SCIP_CALL( SCIPincludePresolTworowbnd(scip) );

   return SCIP_OKAY;
}

/** execution method of presolver */
static
SCIP_DECL_PRESOLEXEC(presolExecTworowbnd)
{  /*lint --e{715}*/
   SCIPMILPMATRIX* matrix;
   SCIP_Bool initialized;
   SCIP_Bool complete;

   assert(result != NULL);
   *result = SCIP_DIDNOTRUN;

   if( (SCIPgetStage(scip) != SCIP_STAGE_PRESOLVING) || SCIPinProbing(scip) || SCIPisNLPEnabled(scip) )
      return SCIP_OKAY;

   if( SCIPgetNVars(scip) == 0 || SCIPisStopped(scip) || SCIPgetNActivePricers(scip) > 0 )
      return SCIP_OKAY;

   if( SCIPgetNContVars(scip)==0 )
      return SCIP_OKAY;

   *result = SCIP_DIDNOTFIND;

   matrix = NULL;
   SCIP_CALL( SCIPmatrixCreate(scip, &matrix, &initialized, &complete) );

   if( initialized && complete )
   {
      int* baserows;
      int nbaserows;
      int ntightenbnds;
      BNDCHGTYPE* tighten;
      int ndeletecons;
      SCIP_Bool* deletecons;
      int ncols;
      int nrows;
      SCIP_Real* lowerbds;
      SCIP_Real* upperbds;

      ncols = SCIPmatrixGetNColumns(matrix);
      nrows = SCIPmatrixGetNRows(matrix);

      ntightenbnds = 0;
      ndeletecons = 0;

      SCIP_CALL( SCIPallocBufferArray(scip, &baserows, nrows) );

      SCIP_CALL( SCIPallocBufferArray(scip, &tighten, ncols) );
      BMSclearMemoryArray(tighten, ncols);

      SCIP_CALL( SCIPallocBufferArray(scip, &lowerbds, ncols) );
      SCIP_CALL( SCIPallocBufferArray(scip, &upperbds, ncols) );
      getBounds(scip, matrix, lowerbds, upperbds);

      SCIP_CALL( SCIPallocBufferArray(scip, &deletecons, nrows) );
      BMSclearMemoryArray(deletecons, nrows);

      SCIP_CALL( getBaseRows(scip, matrix, &nbaserows, baserows) );

      SCIP_CALL( calctworowbnds(scip, matrix,
            nbaserows, baserows, lowerbds, upperbds,
            &ntightenbnds, tighten, &ndeletecons, deletecons) );

      if( ntightenbnds > 0 )
      {
         int c;
         SCIP_VAR* var;
         SCIP_Bool infeas;
         SCIP_Bool tightened;

         for( c = 0; c < ncols; c++ )
         {
            if( tighten[c] == LOWERBOUND )
            {
               var = SCIPmatrixGetVar(matrix, c);
               SCIP_CALL( SCIPtightenVarLb(scip, var, lowerbds[c], FALSE, &infeas, &tightened) );
               if( tightened )
                  ++(*nchgbds);
            }
            else if( tighten[c] == UPPERBOUND )
            {
               var = SCIPmatrixGetVar(matrix, c);
               SCIP_CALL( SCIPtightenVarUb(scip, var, upperbds[c], FALSE, &infeas, &tightened) );
               if( tightened )
                  ++(*nchgbds);
            }
            else if( tighten[c] == BOTHBOUNDS )
            {
               var = SCIPmatrixGetVar(matrix, c);
               SCIP_CALL( SCIPtightenVarLb(scip, var, lowerbds[c], FALSE, &infeas, &tightened) );
               if( tightened )
                  ++(*nchgbds);
               SCIP_CALL( SCIPtightenVarUb(scip, var, upperbds[c], FALSE, &infeas, &tightened) );
               if( tightened )
                  ++(*nchgbds);
            }
         }
      }

      if( ndeletecons > 0 )
      {
         int r;
         for( r = 0; r < nrows; r++ )
         {
            if( deletecons[r] )
            {
               SCIP_CONS* cons;
               cons = SCIPmatrixGetCons(matrix, r);
               SCIP_CALL( SCIPdelCons(scip, cons) );

               (*ndelconss)++;
            }
         }
      }

      SCIPfreeBufferArray(scip, &deletecons);
      SCIPfreeBufferArray(scip, &upperbds);
      SCIPfreeBufferArray(scip, &lowerbds);
      SCIPfreeBufferArray(scip, &tighten);
      SCIPfreeBufferArray(scip, &baserows);
   }

   SCIPmatrixFree(scip, &matrix);

   return SCIP_OKAY;
}

/*
 * presolver specific interface methods
 */

/** creates the tworowbnd presolver and includes it in SCIP */
SCIP_RETCODE SCIPincludePresolTworowbnd(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_PRESOL* presol;

   /* include presolver */
   SCIP_CALL( SCIPincludePresolBasic(scip, &presol, PRESOL_NAME, PRESOL_DESC, PRESOL_PRIORITY, PRESOL_MAXROUNDS,
         PRESOL_TIMING, presolExecTworowbnd, NULL) );
   SCIP_CALL( SCIPsetPresolCopy(scip, presol, presolCopyTworowbnd) );

   return SCIP_OKAY;
}
