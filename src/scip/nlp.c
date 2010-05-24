/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2010 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#pragma ident "@(#) $Id: nlp.c,v 1.2 2010/05/24 17:01:36 bzfviger Exp $"

/**@file   nlp.c
 * @brief  NLP management methods and datastructures
 * @author Thorsten Gellermann
 * @author Stefan Vigerske
 *
 *  In NLP management, we have to differ between the current NLP and the NLPI problem
 *  stored in the NLP solver. All NLP methods affect the current NLP only.
 *  Before solving the current NLP with the NLP solver, the NLP solvers data
 *  has to be updated to the current NLP with a call to nlpFlush().
 *
 *  @todo allow modifications of quadratic and nonquadratic part
 *  @todo allow modifications for rows already in NLP
 *  @todo handle linear rows from LP
 *  @todo replacement of aggregated variables, removal of fixed variables
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/


#include <assert.h>
#include <string.h>

#include "scip/def.h"
#include "scip/message.h"
#include "scip/set.h"
#include "scip/stat.h"
#include "scip/intervalarith.h"
#include "scip/clock.h"
#include "scip/misc.h"
#include "scip/nlp.h"
#include "scip/var.h"
#include "scip/prob.h"
#include "scip/sol.h"
#include "scip/event.h"
#include "nlpi/nlpi.h"
#include "scip/expression.h"
#include "scip/struct_nlp.h"
/* to get value of parameter "nlp/solver" and nlpis array and to get access to set->lp for releasing a variable */
#include "scip/struct_set.h"
/* to get nlp, set, ... in event handling */
#include "scip/struct_scip.h"

/* avoid inclusion of scip.h */
extern
BMS_BLKMEM* SCIPblkmem(
   SCIP*                 scip                /**< SCIP data structure */
   );

/**@name Event Handler declarations */
/**@{ */

#define EVENTHDLR_NAME   "nlpEventHdlr"
#define EVENTHDLR_DESC   "handles all events necessary for maintaining NLP data"

#define ADDNAMESTONLPI   1                   /**< whether to give variable and row names to NLPI */

static
SCIP_DECL_EVENTEXEC( eventExecNlp );
/**@} */

/*
 * private NLP nonlinear row methods
 */

/** announces, that the given coefficient in the constraint matrix changed */
static
void nlrowLinearCoefChanged(
   SCIP_NLROW*           nlrow,              /**< nonlinear row */
   SCIP_VAR*             var,                /**< variables */
   SCIP_NLP*             nlp                 /**< current NLP data */
   )
{
   assert(nlrow != NULL);
   assert(var   != NULL);

   if( nlrow->nlpindex >= -1 )
   {
      /* @todo make this possible */
      SCIPerrorMessage("cannot change coefficients in a row while in NLP\n");
      SCIPABORT();
   }

   nlrow->activity = SCIP_INVALID;
   nlrow->validactivitynlp = -1;
   nlrow->pseudoactivity = SCIP_INVALID;
   nlrow->validpsactivitydomchg = -1;
   nlrow->minactivity = SCIP_INVALID;
   nlrow->maxactivity = SCIP_INVALID;
   nlrow->validactivitybdsdomchg = -1;
}

/** notifies nonlinear row, that its sides were changed */
static
SCIP_RETCODE nlrowSideChanged(
   SCIP_NLROW*           nlrow,              /**< nonlinear row */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_NLP*             nlp,                /**< current NLP data */
   SCIP_SIDETYPE         sidetype            /**< type of side: left or right hand side */
   )
{
   assert(nlrow != NULL);

   if( nlrow->nlpindex >= -1 )
   {
      /* @todo make this possible */
      SCIPerrorMessage("cannot change sides in a row while in NLP\n");
      return SCIP_ERROR;
   }

   return SCIP_OKAY;
}

/** sorts linear part of row entries such that lower variable indices precede higher ones */
static
void nlrowSortLinear(
   SCIP_NLROW*           nlrow                 /**< nonlinear row to be sorted */
   )
{
   assert(nlrow != NULL);

   /* check, if row is already sorted in the LP part, or if the sorting should be delayed */
   if( nlrow->linvarssorted )
      return;

   /* sort linear coefficients */
   SCIPsortPtrRealInt((void**)nlrow->linvars, nlrow->lincoefs, nlrow->linvarsnlpiidx, SCIPvarComp, nlrow->nlinvars);

   nlrow->linvarssorted = TRUE;
}

/** searches linear variable in nonlinear row, returns position in linvars vector or -1 if not found */
static
int nlrowSearchLinearCoef(
   SCIP_NLROW*           nlrow,              /**< nonlinear row to be searched in */
   SCIP_VAR*             var                 /**< variable to be searched for */
   )
{
   int pos;

   assert(nlrow != NULL);
   assert(var   != NULL);

   pos = -1;

   nlrowSortLinear(nlrow);
   SCIPsortedvecFindPtr((void**)nlrow->linvars, SCIPvarComp, (void*)var, nlrow->nlinvars, &pos);

   return pos;
}

/** moves a coefficient in a nonlinear row to a different place, and updates all corresponding data structures */
static
void nlrowMoveLinearCoef(
   SCIP_NLROW*           nlrow,              /**< NLP row */
   int                   oldpos,             /**< old position of coefficient */
   int                   newpos              /**< new position of coefficient */
   )
{
   assert(nlrow != NULL);
   assert(0 <= oldpos && oldpos < nlrow->nlinvars);
   assert(0 <= newpos && newpos < nlrow->nlinvars);
   assert(nlrow->linvars[oldpos] != NULL);

   if( oldpos == newpos )
      return;

   nlrow->linvars[newpos]        = nlrow->linvars[oldpos];
   nlrow->lincoefs[newpos]       = nlrow->lincoefs[oldpos];
   if( nlrow->linvarsnlpiidx != NULL )
      nlrow->linvarsnlpiidx[newpos] = nlrow->linvarsnlpiidx[oldpos];

   /* update sorted flags */
   nlrow->linvarssorted = FALSE;
}

/** adds a previously non existing linear coefficient to a nonlinear row */
static
SCIP_RETCODE nlrowAddLinearCoef(
   SCIP_NLROW*           nlrow,              /**< nonlinear row */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_NLP*             nlp,                /**< current NLP data */
   SCIP_VAR*             var,                /**< variable */
   SCIP_Real             coef                /**< value of coefficient */
   )
{
   int pos;

   assert(nlrow  != NULL);
   assert(blkmem != NULL);
   assert(var    != NULL);
   assert(!SCIPsetIsZero(set, coef));

   SCIP_CALL( SCIPnlrowEnsureLinearSize(nlrow, blkmem, set, nlrow->nlinvars+1) );
   assert(nlrow->linvars  != NULL);
   assert(nlrow->lincoefs != NULL);

   pos = nlrow->nlinvars;
   nlrow->nlinvars++;

   /* insert the variable */
   nlrow->linvars [pos] = var;
   nlrow->lincoefs[pos] = coef;

   nlrowLinearCoefChanged(nlrow, var, nlp);

   /* update sorted flag */
   if( pos > 0 && SCIPvarCompare(nlrow->linvars[pos-1], nlrow->linvars[pos]) > 0 )
      nlrow->linvarssorted = FALSE;

   SCIPdebugMessage("added linear coefficient %g * <%s> at position %d to nonlinear row <%s>\n",
      coef, SCIPvarGetName(var), pos, nlrow->name);

   return SCIP_OKAY;
}

/** deletes coefficient at given position from row */
static
SCIP_RETCODE nlrowDelLinearCoefPos(
   SCIP_NLROW*           nlrow,              /**< nonlinear row to be changed */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_NLP*             nlp,                /**< current NLP data */
   int                   pos                 /**< position in row vector to delete */
   )
{
   SCIP_VAR* var;

   assert(nlrow != NULL);
   assert(set != NULL);
   assert(0 <= pos && pos < nlrow->nlinvars);
   assert(nlrow->linvars[pos] != NULL);

   var = nlrow->linvars[pos];

   /* move last coefficient to position of empty slot (should set sorted flag to FALSE, if not last variable was deleted) */
   nlrowMoveLinearCoef(nlrow, nlrow->nlinvars-1, pos);
   nlrow->nlinvars--;
   assert(pos == nlrow->nlinvars || nlrow->linvarssorted == FALSE);

   nlrowLinearCoefChanged(nlrow, var, nlp);

   return SCIP_OKAY;
}

/** changes a coefficient at given position of a nonlinear row */
static
SCIP_RETCODE nlrowChgLinearCoefPos(
   SCIP_NLROW*           nlrow,              /**< NLP row */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_NLP*             nlp,                /**< current NLP data */
   int                   pos,                /**< position in row vector to change */
   SCIP_Real             coef                /**< new value of coefficient */
   )
{
   assert(nlrow != NULL);
   assert(0 <= pos && pos < nlrow->nlinvars);
   assert(nlrow->linvars[pos] != NULL);

   if( SCIPsetIsZero(set, coef) )
   {
      /* delete existing coefficient */
      SCIP_CALL( nlrowDelLinearCoefPos(nlrow, set, nlp, pos) );
   }
   else if( !SCIPsetIsEQ(set, nlrow->lincoefs[pos], coef) )
   {
      /* change existing coefficient */
      nlrow->lincoefs[pos] = coef;
      nlrowLinearCoefChanged(nlrow, nlrow->linvars[pos], nlp);
   }

   return SCIP_OKAY;
}

#if 0
/** searches quadratic variable in nonlinear row, returns position in quadvars vector or -1 if not found */
static
int nlrowSearchQuadVar(
   SCIP_NLROW*           nlrow,              /**< nonlinear row to be searched in */
   SCIP_VAR*             var                 /**< variable to be searched for */
   )
{
   int i;

   assert(nlrow != NULL);
   assert(var   != NULL);

   /*@todo implement sorting of quadratic coefficients */

   for( i = 0; i < nlrow->nquadvars; ++i )
      if( nlrow->quadvars[i] == var )
         return i;

   return -1;
}
#endif

/** calculates minimal and maximal activity of row w.r.t. the variable's bounds */
static
SCIP_RETCODE nlrowCalcActivityBounds(
   SCIP_NLROW*           nlrow,              /**< nonlinear row */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat                /**< problem statistics data */
   )
{
   SCIP_Real inf;
   SCIP_INTERVAL activity;
   SCIP_INTERVAL bounds;
   int i;

   assert(nlrow != NULL);
   assert(set   != NULL);
   assert(stat  != NULL);

   inf = SCIPsetInfinity(set);

   /* calculate activity bounds */
   SCIPintervalSet(&activity, 0.0);
   for( i = 0; i < nlrow->nlinvars && !SCIPintervalIsEntire(inf, activity); ++i )
   {
      SCIPintervalSetBounds(&bounds, SCIPvarGetLbLocal(nlrow->linvars[i]), SCIPvarGetUbLocal(nlrow->linvars[i]));
      SCIPintervalMulScalar(inf, &bounds, bounds, nlrow->lincoefs[i]);
      SCIPintervalAdd(inf, &activity, activity, bounds);
   }

   /* @todo make sure quadelems is sorted */
   for( i = 0; i < nlrow->nquadelems && !SCIPintervalIsEntire(inf, activity); )
   {
      SCIP_Real a;
      SCIP_INTERVAL b, tmp;
      int idx1;

      idx1 = nlrow->quadelems[i].idx1;
      SCIPintervalSetBounds(&bounds, SCIPvarGetLbLocal(nlrow->quadvars[idx1]), SCIPvarGetUbLocal(nlrow->quadvars[idx1]));

      /* for x_i*(a*x_i + sum_j b_jx_j) we assemble a and sum_j b_jx_j */
      a = 0.0;
      SCIPintervalSet(&b, 0.0);
      do
      {
         if( nlrow->quadelems[i].idx1 == nlrow->quadelems[i].idx2 )
         {
            a = nlrow->quadelems[i].coef;
         }
         else
         {
            SCIPintervalSetBounds(&tmp, SCIPvarGetLbLocal(nlrow->quadvars[nlrow->quadelems[i].idx2]), SCIPvarGetUbLocal(nlrow->quadvars[nlrow->quadelems[i].idx2]));
            SCIPintervalMulScalar(inf, &tmp, tmp, nlrow->quadelems[i].coef);
            SCIPintervalAdd(inf, &b, b, tmp);
         }
         ++i;
      }
      while( i < nlrow->nquadvars && idx1 == nlrow->quadelems[i].idx1 );

      /* compute bounds for a*x_i^2 + b*x_i and add to activity bounds */
      SCIPintervalQuad(inf, &bounds, a, b, bounds);
      SCIPintervalAdd(inf, &activity, activity, bounds);
   }

   if( nlrow->exprtree != NULL && !SCIPintervalIsEntire(inf, activity))
   {
      SCIP_INTERVAL* varvals;
      int n;

      n = SCIPexprtreeGetNVars(nlrow->exprtree);

      SCIP_CALL( SCIPbufferAllocMem(set->buffer, set, (void**)&varvals, n * sizeof(SCIP_INTERVAL)) );

      for( i = 0; i < n; ++i )
      {
         SCIPintervalSetBounds(&varvals[i], SCIPvarGetLbLocal(SCIPexprtreeGetVars(nlrow->exprtree)[i]), SCIPvarGetUbLocal(SCIPexprtreeGetVars(nlrow->exprtree)[i]));
      }

      SCIP_CALL( SCIPexprtreeEvalInt(nlrow->exprtree, inf, varvals, &bounds) );
      SCIPintervalAdd(inf, &activity, activity, bounds);

      SCIPbufferFreeMem(set->buffer, (void**)&varvals, 0);
   }

   nlrow->minactivity = SCIPintervalGetInf(activity);
   nlrow->maxactivity = SCIPintervalGetSup(activity);

   nlrow->validactivitybdsdomchg = stat->domchgcount;

   return SCIP_OKAY;
}

/*
 * public NLP nonlinear row methods
 */

/** create a new nonlinear row
 * the new row is already captured
 */
SCIP_RETCODE SCIPnlrowCreate(
   SCIP_NLROW**          nlrow,              /**< buffer to store pointer to nonlinear row */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   const char*           name,               /**< name of nonlinear row */
   int                   nlinvars,           /**< number of linear variables */
   SCIP_VAR**            linvars,            /**< linear variables, or NULL if nlinvars == 0 */
   SCIP_Real*            lincoefs,           /**< linear coefficients, or NULL if nlinvars == 0 */
   int                   nquadvars,          /**< number of variables in quadratic terms */
   SCIP_VAR**            quadvars,           /**< variables in quadratic terms, or NULL if nquadvars == 0 */
   int                   nquadelems,         /**< number of entries in quadratic term matrix */
   SCIP_QUADELEM*        quadelems,          /**< elements of quadratic term matrix, or NULL if nquadelems == 0 */
   SCIP_EXPRTREE*        exprtree,           /**< expression tree, or NULL */
   SCIP_Real             lhs,                /**< left hand side */
   SCIP_Real             rhs                 /**< right hand side */
   )
{
   assert(nlrow  != NULL);
   assert(blkmem != NULL);
   assert(set    != NULL);
   assert(nlinvars   == 0 || linvars   != NULL);
   assert(nlinvars   == 0 || lincoefs  != NULL);
   assert(nquadvars  == 0 || quadvars  != NULL);
   assert(nquadelems == 0 || quadelems != NULL);
   assert((nquadvars  == 0) == (nquadelems == 0));
   assert(SCIPsetIsLE(set, lhs, rhs));

   SCIP_ALLOC( BMSallocBlockMemory(blkmem, nlrow) );

   /* linear part */
   (*nlrow)->nlinvars = nlinvars;
   (*nlrow)->linvarssize = nlinvars;
   if( nlinvars > 0 )
   {
      SCIP_ALLOC( BMSduplicateBlockMemoryArray(blkmem, &(*nlrow)->linvars,  linvars,  nlinvars) );
      SCIP_ALLOC( BMSduplicateBlockMemoryArray(blkmem, &(*nlrow)->lincoefs, lincoefs, nlinvars) );
      (*nlrow)->linvarssorted = FALSE;
   }
   else
   {
      (*nlrow)->linvars  = NULL;
      (*nlrow)->lincoefs = NULL;
      (*nlrow)->linvarssorted = TRUE;
   }
   (*nlrow)->linvarsnlpiidx = NULL;

   /* quadratic part */
   (*nlrow)->nquadvars  = nquadvars;
   (*nlrow)->nquadelems = nquadelems;
   if( nquadvars > 0 )
   {
      SCIP_ALLOC( BMSduplicateBlockMemoryArray(blkmem, &(*nlrow)->quadvars,  quadvars,  nquadvars) );
      SCIP_ALLOC( BMSduplicateBlockMemoryArray(blkmem, &(*nlrow)->quadelems, quadelems, nquadelems) );
   }
   else
   {
      (*nlrow)->quadvars  = NULL;
      (*nlrow)->quadelems = NULL;
   }
   (*nlrow)->quadvarsnlpiidx = NULL;

   /* nonquadratic part */
   if( exprtree != NULL )
   {
      SCIP_CALL( SCIPexprtreeCopy( blkmem, &(*nlrow)->exprtree, exprtree) );
   }
   else
   {
      (*nlrow)->exprtree = NULL;
   }
   (*nlrow)->exprtreenlpiidx = NULL;

   /* left and right hand sides */
   (*nlrow)->lhs = lhs;
   (*nlrow)->rhs = rhs;

   /* miscellaneous */
   SCIP_ALLOC( BMSduplicateBlockMemoryArray(blkmem, &(*nlrow)->name, name, strlen(name)+1) );
   (*nlrow)->activity = SCIP_INVALID;
   (*nlrow)->validactivitynlp = FALSE;
   (*nlrow)->pseudoactivity = SCIP_INVALID;
   (*nlrow)->validpsactivitydomchg = FALSE;
   (*nlrow)->minactivity = SCIP_INVALID;
   (*nlrow)->maxactivity = SCIP_INVALID;
   (*nlrow)->validactivitybdsdomchg = FALSE;
   (*nlrow)->nlpindex = -2;
   (*nlrow)->nlpiindex = -2;
   (*nlrow)->nuses = 0;

   /* capture the nonlinear row */
   SCIPnlrowCapture(*nlrow);

   return SCIP_OKAY;
}

/** create a nonlinear row that is a copy of a given row
 * the new row is already captured
 */
SCIP_RETCODE SCIPnlrowCreateCopy(
   SCIP_NLROW**          nlrow,              /**< buffer to store pointer to nonlinear row */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_NLROW*           sourcenlrow         /**< nonlinear row to copy */
   )
{
   assert(nlrow  != NULL);
   assert(blkmem != NULL);
   assert(set    != NULL);
   assert(sourcenlrow != NULL);

   SCIP_CALL( SCIPnlrowCreate(nlrow, blkmem, set, sourcenlrow->name,
      sourcenlrow->nlinvars, sourcenlrow->linvars, sourcenlrow->lincoefs,
      sourcenlrow->nquadvars, sourcenlrow->quadvars, sourcenlrow->nquadelems, sourcenlrow->quadelems,
      sourcenlrow->exprtree,
      sourcenlrow->lhs, sourcenlrow->rhs) );

   (*nlrow)->linvarssorted          = sourcenlrow->linvarssorted;
   (*nlrow)->activity               = sourcenlrow->activity;
   (*nlrow)->validactivitynlp       = sourcenlrow->validactivitynlp;
   (*nlrow)->pseudoactivity         = sourcenlrow->pseudoactivity;
   (*nlrow)->validpsactivitydomchg  = sourcenlrow->validpsactivitydomchg;
   (*nlrow)->minactivity            = sourcenlrow->minactivity;
   (*nlrow)->maxactivity            = sourcenlrow->maxactivity;
   (*nlrow)->validactivitybdsdomchg = sourcenlrow->validactivitybdsdomchg;

   return SCIP_OKAY;
}

/** frees a nonlinear row */
SCIP_RETCODE SCIPnlrowFree(
   SCIP_NLROW**          nlrow,              /**< pointer to NLP row */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set                 /**< global SCIP settings */
   )
{
   assert(blkmem != NULL);
   assert(nlrow  != NULL);
   assert(*nlrow != NULL);
   assert((*nlrow)->nuses == 0);
   assert((*nlrow)->nlpindex == -2);
   assert((*nlrow)->nlpiindex == -2);

   /* the NLP should have cleared the following arrays when the row was removed from the NLP */
   assert((*nlrow)->linvarsnlpiidx  == NULL);
   assert((*nlrow)->quadvarsnlpiidx == NULL);
   assert((*nlrow)->exprtreenlpiidx == NULL);

   /* linear part */
   BMSfreeBlockMemoryArrayNull(blkmem, &(*nlrow)->linvars,   (*nlrow)->linvarssize);
   BMSfreeBlockMemoryArrayNull(blkmem, &(*nlrow)->lincoefs,  (*nlrow)->linvarssize);

   /* quadratic part */
   BMSfreeBlockMemoryArrayNull(blkmem, &(*nlrow)->quadvars,  (*nlrow)->nquadvars);
   BMSfreeBlockMemoryArrayNull(blkmem, &(*nlrow)->quadelems, (*nlrow)->nquadelems);

   /* nonquadratic part */
   if( (*nlrow)->exprtree != NULL )
   {
      SCIP_CALL( SCIPexprtreeFree(&(*nlrow)->exprtree) );
   }

   /* miscellaneous */
   BMSfreeBlockMemoryArray(blkmem, &(*nlrow)->name, strlen((*nlrow)->name)+1);

   BMSfreeBlockMemory(blkmem, nlrow);

   return SCIP_OKAY;
}

/** increases usage counter of NLP nonlinear row */
void SCIPnlrowCapture(
   SCIP_NLROW*           nlrow               /**< nonlinear row to capture */
   )
{
   assert(nlrow != NULL);
   assert(nlrow->nuses >= 0);

   SCIPdebugMessage("capture nonlinear row <%s> with nuses=%d\n", nlrow->name, nlrow->nuses);
   nlrow->nuses++;
}

/** decreases usage counter of NLP nonlinear row */
SCIP_RETCODE SCIPnlrowRelease(
   SCIP_NLROW**          nlrow,              /**< nonlinear row to free */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set                 /**< global SCIP settings */
)
{
   assert(blkmem != NULL);
   assert(nlrow  != NULL);
   assert(*nlrow != NULL);
   assert((*nlrow)->nuses >= 1);

   SCIPdebugMessage("release nonlinear row <%s> with nuses=%d\n", (*nlrow)->name, (*nlrow)->nuses);
   (*nlrow)->nuses--;
   if( (*nlrow)->nuses == 0 )
   {
      SCIP_CALL( SCIPnlrowFree(nlrow, blkmem, set) );
   }

   *nlrow = NULL;

   return SCIP_OKAY;
}

/** ensures, that linear coefficient array of nonlinear row can store at least num entries */
SCIP_RETCODE SCIPnlrowEnsureLinearSize(
   SCIP_NLROW*           nlrow,              /**< NLP row */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   int                   num                 /**< minimum number of entries to store */
   )
{
   assert(nlrow != NULL);
   assert(nlrow->nlinvars <= nlrow->linvarssize);

   if( num > nlrow->linvarssize )
   {
      int newsize;

      newsize = SCIPsetCalcMemGrowSize(set, num);
      SCIP_ALLOC( BMSreallocBlockMemoryArray(blkmem, &nlrow->linvars,    nlrow->linvarssize, newsize) );
      SCIP_ALLOC( BMSreallocBlockMemoryArray(blkmem, &nlrow->lincoefs,   nlrow->linvarssize, newsize) );
      if( nlrow->linvarsnlpiidx != NULL )
      {
         SCIP_ALLOC( BMSreallocBlockMemoryArray(blkmem, &nlrow->linvarsnlpiidx, nlrow->linvarssize, newsize) );
      }
      nlrow->linvarssize = newsize;
   }
   assert(num <= nlrow->linvarssize);

   return SCIP_OKAY;
}

/** adds a previously non existing linear coefficient to an NLP nonlinear row */
SCIP_RETCODE SCIPnlrowAddLinearCoef(
   SCIP_NLROW*           nlrow,              /**< NLP nonlinear row */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_NLP*             nlp,                /**< current NLP data */
   SCIP_VAR*             var,                /**< variable */
   SCIP_Real             val                 /**< value of coefficient */
   )
{
   SCIP_CALL( nlrowAddLinearCoef(nlrow, blkmem, set, nlp, var, val) );

   return SCIP_OKAY;
}

/** deletes linear coefficient from nonlinear row */
SCIP_RETCODE SCIPnlrowDelLinearCoef(
   SCIP_NLROW*           nlrow,              /**< nonlinear row to be changed */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_NLP*             nlp,                /**< current NLP data */
   SCIP_VAR*             var                 /**< coefficient to be deleted */
   )
{
   int pos;

   assert(nlrow != NULL);
   assert(var   != NULL);

   /* search the position of the variable in the row's variable vector */
   pos = nlrowSearchLinearCoef(nlrow, var);
   if( pos == -1 )
   {
      SCIPerrorMessage("coefficient for variable <%s> doesn't exist in nonlinear row <%s>\n", SCIPvarGetName(var), nlrow->name);
      return SCIP_INVALIDDATA;
   }
   assert(0 <= pos && pos < nlrow->nlinvars);
   assert(nlrow->linvars[pos] == var);

   /* delete the variable from the row's variable vector */
   SCIP_CALL( nlrowDelLinearCoefPos(nlrow, set, nlp, pos) );

   return SCIP_OKAY;
}

/** changes or adds a linear coefficient to a nonlinear row */
SCIP_RETCODE SCIPnlrowChgLinearCoef(
   SCIP_NLROW*           nlrow,              /**< nonlinear row */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_NLP*             nlp,                /**< current NLP data */
   SCIP_VAR*             var,                /**< variable */
   SCIP_Real             coef                /**< new value of coefficient */
   )
{
   int pos;

   assert(nlrow != NULL);
   assert(nlp != NULL);
   assert(var != NULL);

   /* search the position of the variable in the row's linvars vector */
   pos = nlrowSearchLinearCoef(nlrow, var);

   /* check, if column already exists in the row's linear variables vector */
   if( pos == -1 )
   {
      /* add previously not existing coefficient */
      SCIP_CALL( nlrowAddLinearCoef(nlrow, blkmem, set, nlp, var, coef) );
   }
   else
   {
      /* change the coefficient in the row */
      SCIP_CALL( nlrowChgLinearCoefPos(nlrow, set, nlp, pos, coef) );
   }

   return SCIP_OKAY;
}

/** changes left hand side of nonlinear row */
SCIP_RETCODE SCIPnlrowChgLhs(
   SCIP_NLROW*           nlrow,              /**< nonlinear row */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_NLP*             nlp,                /**< current NLP data */
   SCIP_Real             lhs                 /**< new left hand side */
   )
{
   assert(nlrow != NULL);

   if( !SCIPsetIsEQ(set, nlrow->lhs, lhs) )
   {
      nlrow->lhs = lhs;
      SCIP_CALL( nlrowSideChanged(nlrow, set, nlp, SCIP_SIDETYPE_LEFT) );
   }

   return SCIP_OKAY;
}

/** changes right hand side of nonlinear row */
SCIP_RETCODE SCIPnlrowChgRhs(
   SCIP_NLROW*           nlrow,              /**< nonlinear row */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_NLP*             nlp,                /**< current NLP data */
   SCIP_Real             rhs                 /**< new right hand side */
   )
{
   assert(nlrow != NULL);

   if( !SCIPsetIsEQ(set, nlrow->rhs, rhs) )
   {
      nlrow->rhs = rhs;
      SCIP_CALL( nlrowSideChanged(nlrow, set, nlp, SCIP_SIDETYPE_RIGHT) );
   }

   return SCIP_OKAY;
}

/** recalculates the current activity of a nonlinear row */
SCIP_RETCODE SCIPnlrowRecalcNLPActivity(
   SCIP_NLROW*           nlrow,              /**< nonlinear row */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_NLP*             nlp                 /**< current NLP data */
   )
{
   SCIP_Real val1, val2;
   int i;
   int previdx1;

   assert(nlrow  != NULL);
   assert(stat   != NULL);
   assert(nlp    != NULL);

   if( !SCIPnlpHasSolution(nlp) )
   {
      SCIPerrorMessage("do not have NLP solution for computing NLP activity\n");
      return SCIP_ERROR;
   }
   assert(nlp->primalsolution != NULL);

   nlrow->activity = 0.0;
   for( i = 0; i < nlrow->nlinvars; ++i )
   {
      assert(nlrow->linvars[i] != NULL);

      if( nlrow->linvarsnlpiidx != NULL )
         val1 = nlp->primalsolution[nlp->varmap_nlpi2nlp[nlrow->linvarsnlpiidx[i]]];
      else
         val1 = SCIPnlpGetVarSolVal(nlp, nlrow->linvars[i]);
      assert(val1 < SCIP_INVALID);

      nlrow->activity += nlrow->lincoefs[i] * val1;
   }

   previdx1 = -1;
   for( i = 0; i < nlrow->nquadelems; ++i )
   {
      /* if first index of quadelems is the same as in last round, val1 is still up to date */
      if( previdx1 != nlrow->quadelems[i].idx1 )
      {
         previdx1 = nlrow->quadelems[i].idx1;
         if( nlrow->quadvarsnlpiidx != NULL )
            val1 = nlp->primalsolution[nlp->varmap_nlpi2nlp[nlrow->quadvarsnlpiidx[previdx1]]];
         else
            val1 = SCIPnlpGetVarSolVal(nlp, nlrow->quadvars[previdx1]);
         assert(val1 < SCIP_INVALID);
         if( val1 == 0.0 )
            continue;
      }

      if( nlrow->quadvarsnlpiidx != NULL )
         val2 = nlp->primalsolution[nlp->varmap_nlpi2nlp[nlrow->quadvarsnlpiidx[nlrow->quadelems[i].idx2]]];
      else
         val2 = SCIPnlpGetVarSolVal(nlp, nlrow->quadvars[nlrow->quadelems[i].idx2]);
      assert(val2 < SCIP_INVALID);

      nlrow->activity += nlrow->quadelems[i].coef * val1 * val2;
   }

   if( nlrow->exprtree != NULL )
   {
      SCIP_Real* varvals;
      int n;

      n = SCIPexprtreeGetNVars(nlrow->exprtree);

      SCIP_CALL( SCIPbufferAllocMem(set->buffer, set, (void**)&varvals, n * sizeof(SCIP_Real)) );

      if( nlrow->exprtreenlpiidx != NULL )
         for( i = 0; i < n; ++i )
            varvals[i] = nlp->primalsolution[nlp->varmap_nlpi2nlp[nlrow->exprtreenlpiidx[i]]];
      else
         for( i = 0; i < n; ++i )
            varvals[i] = SCIPnlpGetVarSolVal(nlp, SCIPexprtreeGetVars(nlrow->exprtree)[i]);

      SCIP_CALL( SCIPexprtreeEval(nlrow->exprtree, varvals, &val1) );
      nlrow->activity += val1;

      SCIPbufferFreeMem(set->buffer, (void**)&varvals, 0);
   }

   nlrow->validactivitynlp = stat->nnlps;

   return SCIP_OKAY;
}

/** returns the activity of a nonlinear row in the current NLP solution */
SCIP_RETCODE SCIPnlrowGetNLPActivity(
   SCIP_NLROW*           nlrow,              /**< nonlinear row */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_NLP*             nlp,                /**< current NLP data */
   SCIP_Real*            activity            /**< buffer to store activity value */
   )
{
   assert(nlrow  != NULL);
   assert(stat   != NULL);
   assert(activity != NULL);

   assert(nlrow->validactivitynlp <= stat->nnlps);

   if( nlrow->validactivitynlp != stat->nnlps )
   {
      SCIP_CALL( SCIPnlrowRecalcNLPActivity(nlrow, set, stat, nlp) );
   }
   assert(nlrow->validactivitynlp == stat->nnlps);
   assert(nlrow->activity < SCIP_INVALID);

   *activity = nlrow->activity;

   return SCIP_OKAY;
}

/** gives the feasibility of a nonlinear row in the current NLP solution: negative value means infeasibility */
SCIP_RETCODE SCIPnlrowGetNLPFeasibility(
   SCIP_NLROW*           nlrow,              /**< nonlinear row */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_NLP*             nlp,                /**< current NLP data */
   SCIP_Real*            feasibility         /**< buffer to store feasibility value */
   )
{
   SCIP_Real activity;

   assert(nlrow != NULL);
   assert(feasibility != NULL);

   SCIP_CALL( SCIPnlrowGetNLPActivity(nlrow, set, stat, nlp, &activity) );
   *feasibility = MIN(nlrow->rhs - activity, activity - nlrow->lhs);

   return SCIP_OKAY;
}

/** calculates the current pseudo activity of a nonlinear row */
SCIP_RETCODE SCIPnlrowRecalcPseudoActivity(
   SCIP_NLROW*           nlrow,              /**< nonlinear row */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat                /**< problem statistics */
   )
{
   SCIP_Real val1, val2;
   int i;

   assert(nlrow  != NULL);
   assert(stat   != NULL);

   nlrow->pseudoactivity = 0.0;
   for( i = 0; i < nlrow->nlinvars; ++i )
   {
      assert(nlrow->linvars[i] != NULL);

      val1 = SCIPvarGetBestBound(nlrow->linvars[i]);
      nlrow->pseudoactivity += nlrow->lincoefs[i] * val1;
   }

   for( i = 0; i < nlrow->nquadelems; ++i )
   {
      val1 = SCIPvarGetBestBound(nlrow->quadvars[nlrow->quadelems[i].idx1]);
      if( val1 == 0.0 )
         continue;

      val2 = SCIPvarGetBestBound(nlrow->quadvars[nlrow->quadelems[i].idx2]);
      nlrow->pseudoactivity += nlrow->quadelems[i].coef * val1 * val2;
   }

   if( nlrow->exprtree != NULL )
   {
      SCIP_Real* varvals;
      int n;

      n = SCIPexprtreeGetNVars(nlrow->exprtree);

      SCIP_CALL( SCIPbufferAllocMem(set->buffer, set, (void**)&varvals, n * sizeof(SCIP_Real)) );

      for( i = 0; i < n; ++i )
         varvals[i] = SCIPvarGetBestBound(SCIPexprtreeGetVars(nlrow->exprtree)[i]);

      SCIP_CALL( SCIPexprtreeEval(nlrow->exprtree, varvals, &val1) );
      nlrow->pseudoactivity += val1;

      SCIPbufferFreeMem(set->buffer, (void**)&varvals, 0);
   }

   nlrow->validpsactivitydomchg = stat->domchgcount;

   return SCIP_OKAY;
}

/** returns the pseudo activity of a nonlinear row in the current pseudo solution */
SCIP_RETCODE SCIPnlrowGetPseudoActivity(
   SCIP_NLROW*           nlrow,              /**< nonlinear row */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_Real*            pseudoactivity      /**< buffer to store pseudo activity value */
   )
{
   assert(nlrow != NULL);
   assert(stat  != NULL);
   assert(pseudoactivity != NULL);
   assert(nlrow->validpsactivitydomchg <= stat->domchgcount);

   /* check, if pseudo activity has to be calculated */
   if( nlrow->validpsactivitydomchg != stat->domchgcount )
   {
      SCIP_CALL( SCIPnlrowRecalcPseudoActivity(nlrow, set, stat) );
   }
   assert(nlrow->validpsactivitydomchg == stat->domchgcount);
   assert(nlrow->pseudoactivity < SCIP_INVALID);

   *pseudoactivity = nlrow->pseudoactivity;

   return SCIP_OKAY;
}

/** returns the pseudo feasibility of a nonlinear row in the current pseudo solution: negative value means infeasibility */
SCIP_RETCODE SCIPnlrowGetPseudoFeasibility(
   SCIP_NLROW*           nlrow,              /**< nonlinear row */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_Real*            pseudofeasibility   /**< buffer to store pseudo feasibility value */
   )
{
   SCIP_Real pseudoactivity;

   assert(nlrow != NULL);
   assert(stat  != NULL);
   assert(pseudofeasibility != NULL);

   SCIP_CALL( SCIPnlrowGetPseudoActivity(nlrow, set, stat, &pseudoactivity) );
   *pseudofeasibility = MIN(nlrow->rhs - pseudoactivity, pseudoactivity - nlrow->lhs);

   return SCIP_OKAY;
}

/** returns the activity of a nonlinear row for a given solution */
SCIP_RETCODE SCIPnlrowGetSolActivity(
   SCIP_NLROW*           nlrow,              /**< nonlinear row */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics data */
   SCIP_SOL*             sol,                /**< primal CIP solution */
   SCIP_Real*            activity            /**< buffer to store activity value */
   )
{
   SCIP_Real inf;
   SCIP_Real val1, val2;
   int i;

   assert(nlrow != NULL);
   assert(set != NULL);
   assert(stat != NULL);
   assert(activity != NULL);

   *activity = 0.0;
   for( i = 0; i < nlrow->nlinvars; ++i )
   {
      assert(nlrow->linvars[i] != NULL);

      val1 = SCIPsolGetVal(sol, set, stat, nlrow->linvars[i]);
      if( val1 == SCIP_UNKNOWN ) /*lint !e777*/
      {
         *activity = SCIP_INVALID;
         return SCIP_OKAY;
      }
      *activity += nlrow->lincoefs[i] * val1;
   }

   for( i = 0; i < nlrow->nquadelems; ++i )
   {
      val1 = SCIPsolGetVal(sol, set, stat, nlrow->quadvars[nlrow->quadelems[i].idx1]);
      if( val1 == SCIP_UNKNOWN ) /*lint !e777*/
      {
         *activity = SCIP_INVALID;
         return SCIP_OKAY;
      }
      if( val1 == 0.0 )
         continue;

      val2 = SCIPsolGetVal(sol, set, stat, nlrow->quadvars[nlrow->quadelems[i].idx2]);
      if( val2 == SCIP_UNKNOWN ) /*lint !e777*/
      {
         *activity = SCIP_INVALID;
         return SCIP_OKAY;
      }
      *activity += nlrow->quadelems[i].coef * val1 * val2;
   }

   if( nlrow->exprtree != NULL )
   {
      SCIP_Real* varvals;
      int n;

      n = SCIPexprtreeGetNVars(nlrow->exprtree);

      SCIP_CALL( SCIPbufferAllocMem(set->buffer, set, (void**)&varvals, n * sizeof(SCIP_Real)) );

      for( i = 0; i < n; ++i )
      {
         varvals[i] = SCIPsolGetVal(sol, set, stat, SCIPexprtreeGetVars(nlrow->exprtree)[i]);
         if( val1 == SCIP_UNKNOWN ) /*lint !e777*/
         {
            *activity = SCIP_INVALID;
            SCIPbufferFreeMem(set->buffer, (void**)&varvals, 0);
            return SCIP_OKAY;
         }
      }

      SCIP_CALL( SCIPexprtreeEval(nlrow->exprtree, varvals, &val1) );
      *activity += val1;

      SCIPbufferFreeMem(set->buffer, (void**)&varvals, 0);
   }

   inf = SCIPsetInfinity(set);
   *activity = MAX(*activity, -inf);
   *activity = MIN(*activity, +inf);

   return SCIP_OKAY;
}

/** returns the feasibility of a nonlinear row for the given solution */
SCIP_RETCODE SCIPnlrowGetSolFeasibility(
   SCIP_NLROW*           nlrow,              /**< nonlinear row */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics data */
   SCIP_SOL*             sol,                /**< primal CIP solution */
   SCIP_Real*            feasibility         /**< buffer to store feasibility value */
   )
{
   SCIP_Real activity;

   assert(nlrow != NULL);
   assert(feasibility != NULL);

   SCIP_CALL( SCIPnlrowGetSolActivity(nlrow, set, stat, sol, &activity) );

   *feasibility = MIN(nlrow->rhs - activity, activity - nlrow->lhs);

   return SCIP_OKAY;
}

/** returns the minimal activity of a nonlinear row w.r.t. the variables' bounds */
SCIP_RETCODE SCIPnlrowGetActivityBounds(
   SCIP_NLROW*           nlrow,              /**< nonlinear row */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics data */
   SCIP_Real*            minactivity,        /**< buffer to store minimal activity, or NULL */
   SCIP_Real*            maxactivity         /**< buffer to store maximal activity, or NULL */
   )
{
   assert(nlrow != NULL);
   assert(set != NULL);
   assert(stat != NULL);
   assert(nlrow->validactivitybdsdomchg <= stat->domchgcount);

   /* check, if activity bounds has to be calculated */
   if( nlrow->validactivitybdsdomchg != stat->domchgcount )
   {
      SCIP_CALL( nlrowCalcActivityBounds(nlrow, set, stat) );
   }
   assert(nlrow->validactivitybdsdomchg == stat->domchgcount);
   assert(nlrow->minactivity < SCIP_INVALID);
   assert(nlrow->maxactivity < SCIP_INVALID);

   if( minactivity != NULL )
      *minactivity = nlrow->minactivity;
   if( maxactivity != NULL )
      *maxactivity = nlrow->maxactivity;

   return SCIP_OKAY;
}

/** returns whether the nonlinear row is redundant w.r.t. the variables' bounds */
SCIP_RETCODE SCIPnlrowIsRedundant(
   SCIP_NLROW*           nlrow,              /**< nonlinear row */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics data */
   SCIP_Bool*            isredundant         /**< buffer to store whether row is redundant */
   )
{
   SCIP_Real minactivity;
   SCIP_Real maxactivity;

   assert(nlrow != NULL);
   assert(set != NULL);
   assert(isredundant != NULL);

   SCIP_CALL( SCIPnlrowGetActivityBounds(nlrow, set, stat, &minactivity, &maxactivity) );

   *isredundant = TRUE;
   if( (!SCIPsetIsInfinity(set, -nlrow->lhs) && SCIPsetIsFeasLT(set, minactivity, nlrow->lhs)) ||
       (!SCIPsetIsInfinity(set,  nlrow->rhs) && SCIPsetIsFeasGT(set, maxactivity, nlrow->rhs)) )
      *isredundant = FALSE;

   return SCIP_OKAY;
}

/** output nonlinear row to file stream */
SCIP_RETCODE SCIPnlrowPrint(
   SCIP_NLROW*           nlrow,              /**< NLP row */
   FILE*                 file                /**< output file (or NULL for standard output) */
   )
{
   int i;

   assert(nlrow != NULL);

   /* print row name */
   if( nlrow->name != NULL && nlrow->name[0] != '\0' )
   {
      SCIPmessageFPrintInfo(file, "%s: ", nlrow->name);
   }

   /* print left hand side */
   SCIPmessageFPrintInfo(file, "%.15g <= ", nlrow->lhs);

   /* print coefficients */
   if( nlrow->nlinvars == 0 && nlrow->nquadvars == 0 && nlrow->exprtree == NULL )
      SCIPmessageFPrintInfo(file, "0 ");

   for( i = 0; i < nlrow->nlinvars; ++i )
   {
      assert(nlrow->linvars[i] != NULL);
      assert(SCIPvarGetName(nlrow->linvars[i]) != NULL);
      SCIPmessageFPrintInfo(file, "%+.15g<%s> ", nlrow->lincoefs[i], SCIPvarGetName(nlrow->linvars[i]));
   }

   for( i = 0; i < nlrow->nquadelems; ++i )
   {
      assert(SCIPvarGetName(nlrow->quadvars[nlrow->quadelems[i].idx1]) != NULL);
      assert(SCIPvarGetName(nlrow->quadvars[nlrow->quadelems[i].idx2]) != NULL);
      if( nlrow->quadelems[i].idx1 == nlrow->quadelems[i].idx2 )
         SCIPmessageFPrintInfo(file, "%+.15gsqr(<%s>) ", nlrow->quadelems[i].coef, SCIPvarGetName(nlrow->quadvars[nlrow->quadelems[i].idx1]));
      else
         SCIPmessageFPrintInfo(file, "%+.15g<%s><%s> ", nlrow->quadelems[i].coef, SCIPvarGetName(nlrow->quadvars[nlrow->quadelems[i].idx1]), SCIPvarGetName(nlrow->quadvars[nlrow->quadelems[i].idx2]));
   }

   if( nlrow->exprtree != NULL )
   {
      SCIP_CALL( SCIPexprtreePrintWithNames(nlrow->exprtree, file) );
   }

   /* print right hand side */
   SCIPmessageFPrintInfo(file, "<= %.15g\n", nlrow->rhs);

   return SCIP_OKAY;
}

/** gets number of variables of linear part */
int SCIPnlrowGetNLinearVars(
   SCIP_NLROW*           nlrow               /**< NLP row */
   )
{
   assert(nlrow != NULL);

   return nlrow->nlinvars;
}

/** gets array with variables of linear part */
SCIP_VAR** SCIPnlrowGetLinearVars(
   SCIP_NLROW*           nlrow               /**< NLP row */
   )
{
   assert(nlrow != NULL);

   return nlrow->linvars;
}

/** gets array with coefficients in linear part */
SCIP_Real* SCIPnlrowGetLinearCoefs(
   SCIP_NLROW*           nlrow               /**< NLP row */
   )
{
   assert(nlrow != NULL);

   return nlrow->lincoefs;
}

/** gets number of quadratic variables in quadratic part */
int SCIPnlrowGetNQuadVars(
   SCIP_NLROW*           nlrow               /**< NLP row */
   )
{
   assert(nlrow != NULL);

   return nlrow->nquadvars;
}

/** gets quadratic variables in quadratic part */
SCIP_VAR** SCIPnlrowGetQuadVars(
   SCIP_NLROW*           nlrow               /**< NLP row */
   )
{
   assert(nlrow != NULL);

   return nlrow->quadvars;
}

/** gets number of quadratic elements in quadratic part */
int SCIPnlrowGetNQuadElems(
   SCIP_NLROW*           nlrow               /**< NLP row */
   )
{
   assert(nlrow != NULL);

   return nlrow->nquadelems;
}

/** gets quadratic elements in quadratic part */
SCIP_QUADELEM* SCIPnlrowGetQuadElems(
   SCIP_NLROW*           nlrow               /**< NLP row */
   )
{
   assert(nlrow != NULL);

   return nlrow->quadelems;
}

/** gets array with coefficients in linear part */
void SCIPnlrowGetQuadData(
   SCIP_NLROW*           nlrow,              /**< NLP row */
   int*                  nquadvars,          /**< buffer to store number of variables in quadratic term, or NULL if not of interest */
   SCIP_VAR***           quadvars,           /**< buffer to store pointer to array of variables in quadratic term, or NULL if not of interest */
   int*                  nquadelems,         /**< buffer to store number of entries in quadratic term, or NULL if not of interest */
   SCIP_QUADELEM**       quadelems           /**< buffer to store pointer to arrau of entries in quadratic term, or NULL if not of interest */
   )
{
   assert(nlrow != NULL);

   if( nquadvars  != NULL )
      *nquadvars  = nlrow->nquadvars;
   if( quadvars   != NULL )
      *quadvars   = nlrow->quadvars;
   if( nquadelems != NULL )
      *nquadelems = nlrow->nquadelems;
   if( quadelems  != NULL )
      *quadelems  = nlrow->quadelems;
}

/** gets expression tree */
SCIP_EXPRTREE* SCIPnlrowGetExprtree(
   SCIP_NLROW*           nlrow               /**< NLP row */
   )
{
   assert(nlrow != NULL);

   return nlrow->exprtree;
}

/** returns the left hand side of a nonlinear row */
SCIP_Real SCIPnlrowGetLhs(
   SCIP_NLROW*           nlrow               /**< NLP row */
   )
{
   assert(nlrow != NULL);

   return nlrow->lhs;
}

/** returns the right hand side of a nonlinear row */
SCIP_Real SCIPnlrowGetRhs(
   SCIP_NLROW*           nlrow               /**< NLP row */
   )
{
   assert(nlrow != NULL);

   return nlrow->rhs;
}

/** returns the name of a nonlinear row */
const char* SCIPnlrowGetName(
   SCIP_NLROW*           nlrow               /**< NLP row */
   )
{
   assert(nlrow != NULL);

   return nlrow->name;
}

/** gets position of a nonlinear row in current NLP, or -1 if it is objective, or -2 if not in NLP */
int SCIPnlrowGetNLPPos(
   SCIP_NLROW*           nlrow               /**< NLP row */
   )
{
   assert(nlrow != NULL);

   return nlrow->nlpindex;
}

/** returns TRUE iff row is member of current NLP */
SCIP_Bool SCIPnlrowIsInNLP(
   SCIP_NLROW*           nlrow               /**< NLP row */
   )
{
   assert(nlrow != NULL);

   return nlrow->nlpindex >= -1;
}

/*
 * private NLP methods
 */


/** adds a set of nonlinear rows to the NLP and captures them */
static
SCIP_RETCODE nlpAddNlRows(
   SCIP_NLP*             nlp,                /**< NLP data */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics data */
   int                   nnlrows,            /**< number of nonlinear rows to add */
   SCIP_NLROW**          nlrows              /**< nonlinear rows to add */
   )
{
#ifndef NDEBUG
   int i;
#endif
   int j;
   SCIP_NLROW* nlrow;

   assert(nlp    != NULL);
   assert(blkmem != NULL);
   assert(set    != NULL);
   assert(nlrows != NULL || nnlrows == 0);
   assert(!nlp->indiving);

   SCIP_CALL( SCIPnlpEnsureNlRowsSize(nlp, blkmem, set, nlp->nnlrows + nnlrows) );

   for( j = 0; j < nnlrows; ++j )
   {
      nlrow = nlrows[j];

      /* assert that row is not in NLP (or even NLPI) yet */
      assert(nlrow->nlpindex == -2);
      assert(nlrow->nlpiindex == -2);

#ifndef NDEBUG
      /* assert that variables of row are in NLP */
      for( i = 0; i < nlrow->nlinvars; ++i )
         assert(SCIPhashmapExists(nlp->varhash, nlrow->linvars[i]));

      for( i = 0; i < nlrow->nquadvars; ++i )
         assert(SCIPhashmapExists(nlp->varhash, nlrow->quadvars[i]));

      if( nlrow->exprtree )
      {
         int n;

         n = SCIPexprtreeGetNVars(nlrow->exprtree);
         assert(SCIPexprtreeGetVars(nlrow->exprtree) != NULL || n == 0);

         for( i = 0; i < n; ++i )
            assert(SCIPhashmapExists(nlp->varhash, SCIPexprtreeGetVars(nlrow->exprtree)[i]));
      }
#endif

      /* add row to NLP and capture it */
      nlp->nlrows[nlp->nnlrows + j] = nlrow;
      nlrow->nlpindex = nlp->nnlrows + j;

      SCIPnlrowCapture(nlrow);

      /* if we have a feasible NLP solution and it satisfies the new solution, then it is still feasible
       * if the NLP was globally or locally infeasible, then it stays that way
       * if the NLP was unbounded, then this may not be the case anymore
       */
      if( nlp->solstat <= SCIP_NLPSOLSTAT_FEASIBLE )
      {
         SCIP_Real feasibility;
         SCIP_CALL( SCIPnlrowGetNLPFeasibility(nlrow, set, stat, nlp, &feasibility) );
         if( !SCIPsetIsFeasNegative(set, feasibility) )
            nlp->solstat = SCIP_NLPSOLSTAT_FEASIBLE;
         else
            nlp->solstat = SCIP_NLPSOLSTAT_LOCINFEASIBLE;
      }
      else if( nlp->solstat == SCIP_NLPSOLSTAT_UNBOUNDED )
      {
         nlp->solstat = SCIP_NLPSOLSTAT_UNKNOWN;
      }
   }

   nlp->nnlrows += nnlrows;
   nlp->nunflushednlrowadd += nnlrows;

   return SCIP_OKAY;
}

/** moves a nonlinear row to a different place, and updates all corresponding data structures */
static
void nlpMoveNlrow(
   SCIP_NLP*             nlp,                /**< NLP data structure */
   int                   oldpos,             /**< old position of nonlinear row */
   int                   newpos              /**< new position of nonlinear row */
   )
{
   assert(nlp != NULL);
   assert(0 <= oldpos && oldpos < nlp->nnlrows);
   assert(0 <= newpos && newpos < nlp->nnlrows);
   assert(nlp->nlrows[oldpos] != NULL);

   if( oldpos == newpos )
      return;

   nlp->nlrows[newpos] = nlp->nlrows[oldpos];
   nlp->nlrows[newpos]->nlpindex = newpos;
}

/** deletes nonlinear row with given position from NLP */
static
SCIP_RETCODE nlpDelNlRowPos(
   SCIP_NLP*             nlp,                /**< NLP data structure */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   int                   pos                 /**< position of nonlinear row that is to be removed */
)
{
   SCIP_NLROW* nlrow;

   assert(nlp != NULL);
   assert(blkmem != NULL);
   assert(set != NULL);
   assert(pos >= 0);
   assert(pos < nlp->nnlrows);
   assert(!nlp->indiving);

   nlrow = nlp->nlrows[pos];
   assert(nlrow != NULL);
   assert(nlrow->nlpindex == pos);

   /* clear mapping of variables to NLP indices */
   BMSfreeBlockMemoryArrayNull(blkmem, &nlrow->linvarsnlpiidx,  nlrow->linvarssize);
   BMSfreeBlockMemoryArrayNull(blkmem, &nlrow->quadvarsnlpiidx, nlrow->nquadvars);
   if( nlrow->exprtree != NULL )
      BMSfreeBlockMemoryArrayNull(blkmem, &nlrow->exprtreenlpiidx, SCIPexprtreeGetNVars(nlrow->exprtree));

   /* if row is in NLPI, then mark that it has to be removed in the next flush
    * if row was not in NLPI yet, then we have one unflushed nlrow addition less */
   if( nlrow->nlpiindex >= 0 )
   {
      assert(nlrow->nlpiindex < nlp->nnlrows_solver);
      nlp->nlrowmap_nlpi2nlp[nlrow->nlpiindex] = -1;
      nlrow->nlpiindex = -2;
      ++nlp->nunflushednlrowdel;
   }
   else
   {
      assert(nlrow->nlpiindex == -2); /* if < 0, then -2, since -1 would mean objective function, which makes no sense here */
      --nlp->nunflushednlrowadd;
   }

   /* move NLP row from the end to pos and mark nlrow to be not in NLP anymore */
   nlpMoveNlrow(nlp, nlp->nnlrows-1, pos);
   nlrow->nlpindex = -2;

   /* forget about restriction */
   SCIP_CALL( SCIPnlrowRelease(&nlrow, blkmem, set) );
   --nlp->nnlrows;

   if( nlp->solstat < SCIP_NLPSOLSTAT_LOCOPT )
      nlp->solstat = SCIP_NLPSOLSTAT_FEASIBLE;
   else if( nlp->solstat == SCIP_NLPSOLSTAT_GLOBINFEASIBLE )
      nlp->solstat = SCIP_NLPSOLSTAT_LOCINFEASIBLE;

   return SCIP_OKAY;
}

/** updates bounds on a variable in the NLPI problem */
static
SCIP_RETCODE nlpUpdateVarBounds(
   SCIP_NLP*             nlp,                /**< NLP data */
   SCIP_VAR*             var                 /**< variable which bounds have changed */
   )
{
   int pos;
   SCIP_Real lb;
   SCIP_Real ub;

   assert(nlp != NULL);
   assert(var != NULL);
   assert(nlp->solver != NULL);
   assert(nlp->problem != NULL);
   assert(SCIPhashmapExists(nlp->varhash, var));

   /* original variable bounds are ignored during diving
    * (all variable bounds are reset to their current value in exitDiving) */
   if( nlp->indiving )
      return SCIP_OKAY;

   /* get position of variable in NLP */
   pos = (int) (size_t) SCIPhashmapGetImage(nlp->varhash, var);

   /* if variable not in NLPI yet, nothing to do */
   if( nlp->varmap_nlp2nlpi[pos] == -1 )
      return SCIP_OKAY;

   /* update bounds in NLPI problem */
   pos = nlp->varmap_nlp2nlpi[pos];
   lb = SCIPvarGetLbLocal(var);
   ub = SCIPvarGetUbLocal(var);
   SCIP_CALL( SCIPnlpiChgVarBounds(nlp->solver, nlp->problem, 1, &pos, &lb, &ub) );

   return SCIP_OKAY;
}

/** updates coefficient of a variable in the objective (if its the SCIP objective) */
static
SCIP_RETCODE nlpUpdateScipObjCoef(
   SCIP_NLP*             nlp,                /**< NLP data */
   SCIP_VAR*             var                 /**< variable which bounds have changed */
   )
{
   int pos;
   int objidx;
   SCIP_Real coef;

   assert(nlp != NULL);
   assert(var != NULL);
   assert(nlp->solver != NULL);
   assert(nlp->problem != NULL);
   assert(SCIPhashmapExists(nlp->varhash, var));

   /* if its a user objective, then we have nothing to do here
    * if the objective in the NLPI is not up to date, then we do not need to do something here too */
   if( nlp->objective != NULL || !nlp->objflushed )
      return SCIP_OKAY;

   /* original objective is ignored during diving
    * we just need to remember that at end of diving we have to flush the objective */
   if( nlp->objective == NULL && nlp->indiving )
   {
      nlp->objflushed = FALSE;
      return SCIP_OKAY;
   }

   /* get position of variable in NLP and objective coefficient */
   pos  = (int) (size_t) SCIPhashmapGetImage(nlp->varhash, var);
   coef = SCIPvarGetObj(var);

   /* if variable not in NLPI yet, then we only need to remember to update the objective after variable additions were flushed */
   if( nlp->varmap_nlp2nlpi[pos] == -1 && coef != 0.0 )
   {
      nlp->objflushed = FALSE;
      return SCIP_OKAY;
   }

   /* if we are here, then the objective in the NLPI is up to date,
    * we keep it this way by changing the coefficient of var in the NLPI problem objective */
   pos = nlp->varmap_nlp2nlpi[pos];
   objidx = -1;
   SCIP_CALL( SCIPnlpiChgLinearCoefs(nlp->solver, nlp->problem, objidx, 1, &pos, &coef) );

   return SCIP_OKAY;
}

/** adds new variables to the NLP */
static
SCIP_RETCODE nlpAddVars(
   SCIP_NLP*             nlp,                /**< NLP data structure */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   int                   nvars,              /**< number of variables to add */
   SCIP_VAR**            vars                /**< variable to add to NLP */
   )
{
   int i;
   SCIP_VAR* var;

   assert(nlp    != NULL);
   assert(blkmem != NULL);
   assert(set    != NULL);
   assert(vars   != NULL || nvars == 0);
   assert(!nlp->indiving || nvars == 0);

   if( nvars == 0 )
      return SCIP_OKAY;

   SCIP_CALL( SCIPnlpEnsureVarsSize(nlp, blkmem, set, nlp->nvars + nvars) );
   assert(nlp->sizevars >= nlp->nvars + nvars);

   for( i = 0; i < nvars; ++i )
   {
      var = vars[i];

      assert(SCIPvarIsTransformed(var));
      assert(SCIPvarIsActive(var));
      assert(!SCIPhashmapExists(nlp->varhash, var));

      SCIPvarCapture(var);

      nlp->vars[nlp->nvars+i]            = var;
      nlp->varmap_nlp2nlpi[nlp->nvars+i] = -1;
      SCIP_CALL( SCIPhashmapInsert(nlp->varhash, var, (void*) (size_t) (nlp->nvars+i)) );

      /* update objective, if necessary (new variables have coefficient 0.0 anyway) */
      if( nlp->objective == NULL && SCIPvarGetObj(var) != 0.0 )
      {
         SCIP_CALL( nlpUpdateScipObjCoef(nlp, var) );
      }

      /* let's keep the previous initial guess and set it for the new variable to the best bound (if SCIP objective) or 0.0 projected on bounds (if user objective)
       * (since there can be no row that uses this variable yet, this seems a good guess) */
      if( nlp->haveinitguess )
      {
         assert(nlp->initialguess != NULL);

         if( nlp->objective == NULL )
            nlp->initialguess[nlp->nvars+i] = SCIPvarGetBestBound(var);
         else
            nlp->initialguess[nlp->nvars+i] = MIN(SCIPvarGetUbLocal(var), MAX(SCIPvarGetLbLocal(var), 0.0));
      }

      /* if we have a feasible NLP solution, then it remains feasible
       * but if we use the SCIP objective, then we have to update the objective function
       */
      if( nlp->solstat <= SCIP_NLPSOLSTAT_FEASIBLE )
      {
         assert(nlp->primalsolution != NULL);

         if( nlp->objective == NULL )
         {
            nlp->primalsolution[nlp->nvars+i] = SCIPvarGetBestBound(var);
            nlp->primalsolobjval += SCIPvarGetObj(var) * nlp->primalsolution[nlp->nvars+i];
         }
         else
         {
            nlp->primalsolution[nlp->nvars+i] = 0.0;
         }
         nlp->solstat = SCIP_NLPSOLSTAT_FEASIBLE;
      }

      /* catch events on variable */
      SCIP_CALL( SCIPvarCatchEvent(var, blkmem, set,
         SCIP_EVENTTYPE_VARFIXED | SCIP_EVENTTYPE_BOUNDCHANGED | SCIP_EVENTTYPE_OBJCHANGED,
         nlp->eventhdlr, (SCIP_EVENTDATA*)nlp, NULL) ); /* @todo should store event filter position is nlp? */
   }

   nlp->nvars += nvars;
   nlp->nunflushedvaradd += nvars;

   return SCIP_OKAY;
}

/** moves a variable to a different place, and updates all corresponding data structures */
static
SCIP_RETCODE nlpMoveVar(
   SCIP_NLP*             nlp,                /**< NLP data structure */
   int                   oldpos,             /**< old position of variable */
   int                   newpos              /**< new position of variable */
   )
{
   int nlpipos;

   assert(nlp != NULL);
   assert(0 <= oldpos && oldpos < nlp->nvars);
   assert(0 <= newpos && newpos < nlp->nvars);
   assert(nlp->vars[oldpos] != NULL);

   if( oldpos == newpos )
      return SCIP_OKAY;

   SCIP_CALL( SCIPhashmapSetImage(nlp->varhash, nlp->vars[oldpos], (void*) (size_t) newpos) );
   nlp->vars[newpos]            = nlp->vars[oldpos];
   nlp->varmap_nlp2nlpi[newpos] = nlp->varmap_nlp2nlpi[oldpos];
   if( nlp->primalsolution != NULL )
      nlp->primalsolution[newpos] = nlp->primalsolution[oldpos];
   if( nlp->initialguess != NULL )
      nlp->initialguess[newpos] = nlp->initialguess[oldpos];

   nlpipos = nlp->varmap_nlp2nlpi[newpos];
   if( nlpipos > 0 )
      nlp->varmap_nlpi2nlp[nlpipos] = newpos;

   return SCIP_OKAY;
}

/** deletes variable with given position from NLP */
static
SCIP_RETCODE nlpDelVarPos(
   SCIP_NLP*             nlp,                /**< NLP data structure */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_LP*              lp,                 /**< SCIP LP, needed if a column-variable is freed */
   int                   pos                 /**< position of nonlinear row that is to be removed */
)
{
   SCIP_VAR* var;
#ifndef NDEBUG
   int i;
#endif
   int nlpipos;

   assert(nlp    != NULL);
   assert(blkmem != NULL);
   assert(set    != NULL);
   assert(pos >= 0);
   assert(pos < nlp->nvars);
   assert(!nlp->indiving);

   var = nlp->vars[pos];
   assert(var != NULL);

#ifndef NDEBUG
   /* assert that variable is not used by any nonlinear row */
   for( i = 0; i < nlp->nnlrows; ++i )
   {
      int j;
      SCIP_NLROW* nlrow;

      nlrow = nlp->nlrows[i];
      assert(nlrow != NULL);

      /* use nlrowSearchLinearCoef only if already sorted, since otherwise we may change the solving process slightly */
      if( nlrow->linvarssorted )
         assert( nlrowSearchLinearCoef(nlrow, var) == -1 );
      else
         for( j = 0; j < nlrow->nlinvars; ++j )
            assert( nlrow->linvars[j] != var );

      for( j = 0; j < nlrow->nquadvars; ++j )
         assert( nlrow->quadvars[j] != var );

      assert(nlrow->exprtree == NULL || SCIPexprtreeFindVar(nlrow->exprtree, var) == -1);
   }
#endif

   /* if we had a feasible solution and used the SCIP objective, then adjust objective function value
    * if NLP was unbounded before, then maybe it is not anymore */
   if( nlp->solstat <= SCIP_NLPSOLSTAT_FEASIBLE )
   {
      assert(nlp->primalsolution != NULL);
      if( nlp->objective == NULL )
         nlp->primalsolobjval -= SCIPvarGetObj(var) * nlp->primalsolution[pos];
   }
   else if( nlp->solstat == SCIP_NLPSOLSTAT_UNBOUNDED )
      nlp->solstat = SCIP_NLPSOLSTAT_UNKNOWN;

   /* if variable is in NLPI problem, mark that we have to remember to delete it there
    * if it was not in the NLPI yet, then we have one unflushed var addition less now */
   nlpipos = nlp->varmap_nlp2nlpi[pos];
   if( nlpipos >= 0 )
   {
      assert(nlpipos < nlp->nvars_solver);

      nlp->varmap_nlpi2nlp[nlpipos] = -1;
      ++nlp->nunflushedvardel;
   }
   else
      --nlp->nunflushedvaradd;

   /* drop events on variable */
   SCIP_CALL( SCIPvarDropEvent(var, blkmem, set,
      SCIP_EVENTTYPE_VARFIXED | SCIP_EVENTTYPE_BOUNDCHANGED | SCIP_EVENTTYPE_OBJCHANGED,
      nlp->eventhdlr, (SCIP_EVENTDATA*)nlp, -1) );

   /* move variable from end to pos */
   SCIP_CALL( nlpMoveVar(nlp, nlp->nvars-1, pos) );

   /* forget about variable */
   SCIP_CALL( SCIPhashmapRemove(nlp->varhash, var) );
   SCIP_CALL( SCIPvarRelease(&var, blkmem, set, lp) );
   --nlp->nvars;

   return SCIP_OKAY;
}

/** gets NLPI variable indices of variable in nonlinear row */
static
SCIP_RETCODE nlpSetupNlpiIndices(
   SCIP_NLP*             nlp,                /**< NLP data */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_NLROW*           nlrow               /**< nonlinear row */
   )
{
   int i;
   SCIP_VAR* var;

   assert(nlp    != NULL);
   assert(blkmem != NULL);
   assert(set    != NULL);
   assert(nlrow  != NULL);
   assert(nlrow->nlpiindex == -2); /* row should not be in NLPI yet */
   assert(nlrow->linvarsnlpiidx  == NULL);
   assert(nlrow->quadvarsnlpiidx == NULL);
   assert(nlrow->exprtreenlpiidx == NULL);

   /* get indices of variables in linear part of row */
   if( nlrow->nlinvars > 0 )
   {
      assert(nlrow->linvars  != NULL);
      assert(nlrow->lincoefs != NULL);

      if( nlrow->linvarsnlpiidx == NULL )
      {
         SCIP_ALLOC( BMSallocBlockMemoryArray(blkmem, &nlrow->linvarsnlpiidx, nlrow->linvarssize) );
      }

      for( i = 0; i < nlrow->nlinvars; ++i )
      {
         var = nlrow->linvars[i];
         assert(var != NULL);

         assert(SCIPhashmapExists(nlp->varhash, var));
         nlrow->linvarsnlpiidx[i] = (size_t) (void*) SCIPhashmapGetImage(nlp->varhash, var);
      }
   }

   /* get indices of variables in quadratic part of row */
   if( nlrow->nquadvars > 0 )
   {
      assert(nlrow->quadvars    != NULL);
      assert(nlrow->nquadelems  > 0);
      assert(nlrow->quadelems   != NULL);

      if( nlrow->quadvarsnlpiidx == NULL )
      {
         SCIP_ALLOC( BMSallocBlockMemoryArray(blkmem, &nlrow->quadvarsnlpiidx, nlrow->nquadvars) );
      }

      for( i = 0; i < nlrow->nquadvars; ++i )
      {
         var = nlrow->quadvars[i];
         assert(var != NULL);

         assert(SCIPhashmapExists(nlp->varhash, var));
         nlrow->quadvarsnlpiidx[i] = (size_t) (void*) SCIPhashmapGetImage(nlp->varhash, var);
      }
   }

   /* get indices of variables in expression tree part of row */
   if( nlrow->exprtree != NULL )
   {
      int n;

      n = SCIPexprtreeGetNVars(nlrow->exprtree);
      assert(n == 0 || SCIPexprtreeGetVars(nlrow->exprtree) != NULL);

      if( nlrow->exprtreenlpiidx == NULL )
      {
         SCIP_ALLOC( BMSallocBlockMemoryArray(blkmem, &nlrow->exprtreenlpiidx, n) );
      }

      for( i = 0; i < n; ++i )
      {
         var = SCIPexprtreeGetVars(nlrow->exprtree)[i];
         assert(var != NULL);

         assert(SCIPhashmapExists(nlp->varhash, var));
         nlrow->exprtreenlpiidx[i] = (size_t) (void*) SCIPhashmapGetImage(nlp->varhash, var);
      }
   }

   return SCIP_OKAY;
}

/** ensures, that NLPI variables array of NLP can store at least num entries */
static
SCIP_RETCODE nlpEnsureVarsSolverSize(
   SCIP_NLP*             nlp,                /**< NLP data */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   int                   num                 /**< minimum number of entries to store */
   )
{
   assert(nlp    != NULL);
   assert(blkmem != NULL);
   assert(set    != NULL);
   assert(nlp->nvars_solver <= nlp->sizevars_solver);

   if( num > nlp->sizevars_solver )
   {
      int newsize;

      newsize = SCIPsetCalcMemGrowSize(set, num);
      SCIP_ALLOC( BMSreallocBlockMemoryArray(blkmem, &nlp->varmap_nlpi2nlp, nlp->sizevars_solver, newsize) );

      nlp->sizevars_solver = newsize;
   }
   assert(num <= nlp->sizevars_solver);

   return SCIP_OKAY;
}

/** ensures, that NLPI nonlinear rows array of NLP can store at least num entries */
static
SCIP_RETCODE nlpEnsureNlRowsSolverSize(
   SCIP_NLP*             nlp,                /**< NLP data */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   int                   num                 /**< minimum number of entries to store */
   )
{
   assert(nlp    != NULL);
   assert(blkmem != NULL);
   assert(set    != NULL);
   assert(nlp->nnlrows_solver <= nlp->sizenlrows_solver);

   if( num > nlp->sizenlrows_solver )
   {
      int newsize;

      newsize = SCIPsetCalcMemGrowSize(set, num);
      SCIP_ALLOC( BMSreallocBlockMemoryArray(blkmem, &nlp->nlrowmap_nlpi2nlp, nlp->sizenlrows_solver, newsize) );

      nlp->sizenlrows_solver = newsize;
   }
   assert(num <= nlp->sizenlrows_solver);

   return SCIP_OKAY;
}

/** deletes rows from the NLPI problem that have been marked as to remove */
static
SCIP_RETCODE nlpFlushNlRowDeletions(
   SCIP_NLP*             nlp,                /**< NLP data */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set                 /**< global SCIP settings */
)
{
   int         j;
   int         c;      /* counts the number of rows to delete */
   int*        rowset; /* marks which rows to delete and stores new indices */
   SCIP_NLROW* nlrow;

   assert(nlp    != NULL);
   assert(blkmem != NULL);
   assert(set    != NULL);
   assert(nlp->solver != NULL);
   assert(nlp->problem != NULL);
   assert(nlp->nunflushednlrowdel >= 0);
   assert(!nlp->indiving);

   if( nlp->nunflushednlrowdel == 0 )
   {
#ifndef NDEBUG
      /* check that there are really no pending removals of nonlinear rows */
      for( j = 0; j < nlp->nnlrows_solver; ++j )
         assert(nlp->nlrowmap_nlpi2nlp[j] >= 0);
#endif
      return SCIP_OKAY;
   }

   /* create marker which rows have to be deleted */
   SCIP_CALL( SCIPbufferAllocMem(set->buffer, set, (void**)&rowset, nlp->nnlrows_solver * sizeof(int)) );
   c = 0;
   for( j = 0; j < nlp->nnlrows_solver; ++j )
   {
      if( nlp->nlrowmap_nlpi2nlp[j] == -1 )
      {
         rowset[j] = 1;
         ++c;
      }
      else
         rowset[j] = 0;
   }
   assert(c == nlp->nunflushednlrowdel);

   /* remove rows from NLPI problem */
   SCIP_CALL( SCIPnlpiDelConsSet(nlp->solver, nlp->problem, rowset) );

   /* update NLPI row indices */
   for( j = 0; j < nlp->nnlrows_solver; ++j )
   {
      assert(rowset[j] <= j); /* we assume that the NLP solver did not move a row behind its previous position!! */
      if( rowset[j] < 0 )
      {
         /* assert that row was marked as deleted */
         assert(nlp->nlrowmap_nlpi2nlp[j] == -1);
      }
      else if( rowset[j] < j )
      {
         /* nlrow at position j moved (forward) to position rowset[j] */
         assert(nlp->nlrowmap_nlpi2nlp[j] >= 0);
         assert(nlp->nlrowmap_nlpi2nlp[j] < nlp->nnlrows);

         nlrow = nlp->nlrows[nlp->nlrowmap_nlpi2nlp[j]];
         assert(nlrow->nlpiindex == j);

         /* there should be no row at the new position already */
         assert(nlp->nlrowmap_nlpi2nlp[rowset[j]] == -1);

         nlrow->nlpiindex = rowset[j];
         nlp->nlrowmap_nlpi2nlp[rowset[j]] = nlrow->nlpindex;
      }
      else
      {
         /* row j stays at position j */
         assert(nlp->nlrowmap_nlpi2nlp[j] >= 0);
         assert(nlp->nlrowmap_nlpi2nlp[j] < nlp->nnlrows);
         assert(nlp->nlrows[nlp->nlrowmap_nlpi2nlp[j]]->nlpiindex == j);
      }
   }
   nlp->nnlrows_solver -= c;
   nlp->nunflushednlrowdel = 0;

   /* cleanup */
   SCIPbufferFreeMem(set->buffer, (void**)&rowset, 0);

   return SCIP_OKAY;
}

/** deletes variables from the NLPI problem that have been marked as to remove
 * assumes that there are no pending row deletions (nlpFlushNlRowDeletions should be called first)
 */
static
SCIP_RETCODE nlpFlushVarDeletions(
   SCIP_NLP*             nlp,                /**< NLP data */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set                 /**< global SCIP settings */
)
{
   int  i;
   int  c;      /* counter on number of variables to remove in solver */
   int* colset; /* marks which variables to delete and stores new indices */

   assert(nlp    != NULL);
   assert(blkmem != NULL);
   assert(set    != NULL);
   assert(nlp->solver != NULL);
   assert(nlp->problem != NULL);
   assert(nlp->nunflushedvardel >= 0);
   assert(nlp->nunflushednlrowdel == 0);
   assert(!nlp->indiving);

   if( nlp->nunflushedvardel == 0 )
   {
#ifndef NDEBUG
      /* check that there are really no pending removals of variables */
      for( i = 0; i < nlp->nvars_solver; ++i )
         assert(nlp->varmap_nlpi2nlp[i] >= 0);
#endif
      return SCIP_OKAY;
   }

   /* create marker which variables have to be deleted */
   SCIP_CALL( SCIPbufferAllocMem(set->buffer, set, (void**)&colset, nlp->nvars_solver * sizeof(int)) );
   c = 0;
   for( i = 0; i < nlp->nvars_solver; ++i )
   {
      if( nlp->varmap_nlpi2nlp[i] == -1 )
      {
         colset[i] = 1;
         ++c;
      }
      else
         colset[i] = 0;
   }
   assert(c == nlp->nunflushedvardel);

   /* delete variables from NLPI problem */
   SCIP_CALL( SCIPnlpiDelVarSet(nlp->solver, nlp->problem, colset) );

   /* update NLPI variable indices */
   for( i = 0; i < nlp->nvars_solver; ++i )
   {
      assert(colset[i] <= i); /* we assume that the NLP solver did not move a variable behind its previous position!! */
      if( colset[i] < 0 )
      {
         /* assert that variable was marked as deleted */
         assert(nlp->varmap_nlpi2nlp[i] == -1);
      }
      else if( colset[i] < i)
      {
         /* variable at position i moved (forward) to position colset[i] */
         int varpos;

         varpos = nlp->varmap_nlpi2nlp[i]; /* position of variable i in NLP */
         assert(varpos >= 0);
         assert(varpos < nlp->nvars);
         assert(nlp->varmap_nlp2nlpi[varpos] == i);

         /* there should be no variable at the new position already */
         assert(nlp->varmap_nlpi2nlp[colset[i]] == -1);

         nlp->varmap_nlp2nlpi[varpos] = colset[i];
         nlp->varmap_nlpi2nlp[colset[i]] = varpos;
      }
      else
      {
         /* variable i stays at position i */
         assert(nlp->varmap_nlpi2nlp[i] >= 0);
         assert(nlp->varmap_nlpi2nlp[i] < nlp->nvars);
         assert(nlp->varmap_nlp2nlpi[nlp->varmap_nlpi2nlp[i]] == i);
      }
   }

   nlp->nvars_solver -= c;
   nlp->nunflushedvardel = 0;

   /* cleanup */
   SCIPbufferFreeMem(set->buffer, (void**)&colset, 0);

   return SCIP_OKAY;
}

/** adds nonlinear rows to NLPI problem that have been added to NLP before
 * assumes that there are no pending variable additions or deletions (nlpFlushVarDeletions and nlpFlushVarAdditions should be called first) */
static
SCIP_RETCODE nlpFlushNlRowAdditions(
   SCIP_NLP*             nlp,                /**< NLP data */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set                 /**< global SCIP settings */
)
{
   int c, i, j;
   SCIP_Bool   havequad;
   SCIP_NLROW* nlrow;
   SCIP_Real*  lhss;
   SCIP_Real*  rhss;
   int*        nlinvars;
   int**       linidxs;
   SCIP_Real** lincoefs;
   int*        nquadelems;
   SCIP_QUADELEM** quadelems;
   int**       nlidxs;
   SCIP_EXPRTREE** exprtrees;
   const char** names;

   assert(nlp    != NULL);
   assert(blkmem != NULL);
   assert(set    != NULL);
   assert(nlp->nunflushednlrowadd >= 0);
   assert(nlp->nunflushedvaradd == 0);
   assert(nlp->nunflushedvardel == 0);
   assert(!nlp->indiving);

   if( nlp->nunflushednlrowadd == 0 )
   {
#ifndef NDEBUG
      /* check that there are really no pending additions of variables */
      for( i = 0; i < nlp->nnlrows; ++i )
         assert(nlp->nlrows[i]->nlpiindex >= 0);
#endif
      return SCIP_OKAY;
   }

   SCIP_CALL( nlpEnsureNlRowsSolverSize(nlp, blkmem, set, nlp->nnlrows_solver + nlp->nunflushednlrowadd) );

   SCIP_CALL( SCIPbufferAllocMem(set->buffer, set, (void**)&lhss,        nlp->nunflushednlrowadd * sizeof(SCIP_Real)) );
   SCIP_CALL( SCIPbufferAllocMem(set->buffer, set, (void**)&rhss,        nlp->nunflushednlrowadd * sizeof(SCIP_Real)) );
   SCIP_CALL( SCIPbufferAllocMem(set->buffer, set, (void**)&nlinvars,    nlp->nunflushednlrowadd * sizeof(int)) );
   SCIP_CALL( SCIPbufferAllocMem(set->buffer, set, (void**)&linidxs,     nlp->nunflushednlrowadd * sizeof(int*)) );
   SCIP_CALL( SCIPbufferAllocMem(set->buffer, set, (void**)&lincoefs,    nlp->nunflushednlrowadd * sizeof(SCIP_Real*)) );
   SCIP_CALL( SCIPbufferAllocMem(set->buffer, set, (void**)&nquadelems,  nlp->nunflushednlrowadd * sizeof(int)) );
   SCIP_CALL( SCIPbufferAllocMem(set->buffer, set, (void**)&quadelems,   nlp->nunflushednlrowadd * sizeof(SCIP_QUADELEM*)) );
   SCIP_CALL( SCIPbufferAllocMem(set->buffer, set, (void**)&nlidxs,      nlp->nunflushednlrowadd * sizeof(int*)) );
   SCIP_CALL( SCIPbufferAllocMem(set->buffer, set, (void**)&exprtrees,   nlp->nunflushednlrowadd * sizeof(SCIP_EXPRTREE*)) );
#if ADDNAMESTONLPI
   SCIP_CALL( SCIPbufferAllocMem(set->buffer, set, (void**)&names,       nlp->nunflushednlrowadd * sizeof(const char*)) );
#else
   names = NULL;
#endif

   c = 0;
   havequad = FALSE; /* indicates whether some row with a nonzero quadratic part was added */
   for( i = 0; i < nlp->nnlrows; ++i )
   {
      nlrow = nlp->nlrows[i];
      assert(nlrow != NULL);

      /* skip nonlinear rows already in NLPI problem */
      if( nlrow->nlpiindex >= 0 )
         continue;
      assert(c < nlp->nunflushednlrowadd);

      SCIP_CALL( nlpSetupNlpiIndices(nlp, blkmem, set, nlrow) );
      assert(nlrow->linvarsnlpiidx  != NULL || nlrow->nlinvars  == 0);
      assert(nlrow->quadvarsnlpiidx != NULL || nlrow->nquadvars == 0);
      assert(nlrow->exprtreenlpiidx != NULL || nlrow->exprtree  == NULL);

      nlp->nlrowmap_nlpi2nlp[nlp->nnlrows_solver+c] = i;
      nlrow->nlpiindex = nlp->nnlrows_solver+c;

      lhss[c]   = nlrow->lhs;
      rhss[c]   = nlrow->rhs;

      nlinvars[c]    = nlrow->nlinvars;
      linidxs[c]     = nlrow->linvarsnlpiidx;
      lincoefs[c]    = nlrow->lincoefs;

      nquadelems[c]  = nlrow->nquadelems;
      if( nlrow->nquadelems > 0 )
      {
         SCIP_CALL( SCIPbufferAllocMem(set->buffer, set, (void**)&quadelems[c], nlrow->nquadelems * sizeof(SCIP_QUADELEM)) );
         for( j = 0; j < nlrow->nquadelems; ++j )
         {
            quadelems[c][j].idx1 = MIN(nlrow->quadvarsnlpiidx[nlrow->quadelems[j].idx1], nlrow->quadvarsnlpiidx[nlrow->quadelems[j].idx2]);
            quadelems[c][j].idx2 = MAX(nlrow->quadvarsnlpiidx[nlrow->quadelems[j].idx1], nlrow->quadvarsnlpiidx[nlrow->quadelems[j].idx2]);
            quadelems[c][j].coef = nlrow->quadelems[j].coef;
         }
         havequad = TRUE;
      }
      else
         quadelems[c] = NULL;

      nlidxs[c]      = nlrow->exprtreenlpiidx;
      exprtrees[c]   = nlrow->exprtree;

#if ADDNAMESTONLPI
      names[c]       = nlrow->name;
#endif

      ++c;

#ifdef NDEBUG
      /* have c vars to add already, there can be no more */
      if( c == nlp->nunflushednlrowadd )
         break;
#endif
   }
   assert(c == nlp->nunflushednlrowadd);

   nlp->nnlrows_solver += c;

   SCIP_CALL( SCIPnlpiAddConstraints(nlp->solver, nlp->problem, c, lhss, rhss,
      nlinvars, linidxs, lincoefs,
      nquadelems, quadelems,
      nlidxs, exprtrees,
      names) );

#if ADDNAMESTONLPI
   SCIPbufferFreeMem(set->buffer, (void**)&names, 0);
#endif
   SCIPbufferFreeMem(set->buffer, (void**)&lhss, 0);
   SCIPbufferFreeMem(set->buffer, (void**)&rhss, 0);
   SCIPbufferFreeMem(set->buffer, (void**)&nlinvars, 0);
   SCIPbufferFreeMem(set->buffer, (void**)&linidxs, 0);
   SCIPbufferFreeMem(set->buffer, (void**)&lincoefs, 0);
   if( havequad )
      for( c = 0; c < nlp->nunflushednlrowadd; ++c )
         if( quadelems[c] != NULL )
            SCIPbufferFreeMem(set->buffer, (void**)&quadelems[c], 0);
   SCIPbufferFreeMem(set->buffer, (void**)&nquadelems, 0);
   SCIPbufferFreeMem(set->buffer, (void**)&quadelems, 0);
   SCIPbufferFreeMem(set->buffer, (void**)&nlidxs, 0);
   SCIPbufferFreeMem(set->buffer, (void**)&exprtrees, 0);

   nlp->nunflushednlrowadd = 0;

   return SCIP_OKAY;
}


/** adds variables to NLPI problem that have been added to NLP before
 * may set nlp->objflushed to TRUE if objective is SCIP objective and a variable with nonzero obj.coefficient is added to the NLPI problem */
static
SCIP_RETCODE nlpFlushVarAdditions(
   SCIP_NLP*             nlp,                /**< NLP data */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set                 /**< global SCIP settings */
)
{
   int i, c;
   SCIP_Real*  lbs;
   SCIP_Real*  ubs;
   const char** names;

   assert(nlp    != NULL);
   assert(blkmem != NULL);
   assert(set    != NULL);
   assert(nlp->solver != NULL);
   assert(nlp->problem != NULL);
   assert(nlp->nunflushedvaradd >= 0);
   assert(!nlp->indiving);

   if( nlp->nunflushedvaradd == 0 )
   {
#ifndef NDEBUG
      /* check that there are really no pending additions of variables */
      for( i = 0; i < nlp->nvars; ++i )
         assert(nlp->varmap_nlp2nlpi[i] >= 0);
#endif
      return SCIP_OKAY;
   }

   SCIP_CALL( nlpEnsureVarsSolverSize(nlp, blkmem, set, nlp->nvars_solver + nlp->nunflushedvaradd) );

   SCIP_CALL( SCIPbufferAllocMem(set->buffer, set, (void**)&lbs,   nlp->nunflushedvaradd * sizeof(SCIP_Real)) );
   SCIP_CALL( SCIPbufferAllocMem(set->buffer, set, (void**)&ubs,   nlp->nunflushedvaradd * sizeof(SCIP_Real)) );
#if ADDNAMESTONLPI
   SCIP_CALL( SCIPbufferAllocMem(set->buffer, set, (void**)&names, nlp->nunflushedvaradd * sizeof(const char*)) );
#else
   names = NULL;
#endif

   c = 0;
   for( i = 0; i < nlp->nvars; ++i )
   {
      /* skip variables already in NLPI problem */
      if( nlp->varmap_nlp2nlpi[i] >= 0 )
         continue;
      assert(c < nlp->nunflushedvaradd);

      nlp->varmap_nlpi2nlp[nlp->nvars_solver+c] = i;
      nlp->varmap_nlp2nlpi[i] = nlp->nvars_solver+c;
      lbs[c]   = SCIPvarGetLbLocal(nlp->vars[i]);
      ubs[c]   = SCIPvarGetUbLocal(nlp->vars[i]);
#if ADDNAMESTONLPI
      names[c] = SCIPvarGetName(nlp->vars[i]);
#endif
      ++c;

      /* if we use the SCIP objective function and the new variable has a nonzero objective coefficient,
       * then the objective need to be updated */
      if( nlp->objective == NULL && !SCIPsetIsZero(set, SCIPvarGetObj(nlp->vars[i])) )
         nlp->objflushed = FALSE;

#ifdef NDEBUG
      /* have c vars to add already, there can be no more */
      if( c == nlp->nunflushedvaradd )
         break;
#endif
   }
   assert(c == nlp->nunflushedvaradd);

   nlp->nvars_solver += c;

   SCIP_CALL( SCIPnlpiAddVars(nlp->solver, nlp->problem, c, lbs, ubs, names) );

#if ADDNAMESTONLPI
   SCIPbufferFreeMem(set->buffer, (void**)&names, 0);
#endif
   SCIPbufferFreeMem(set->buffer, (void**)&lbs, 0);
   SCIPbufferFreeMem(set->buffer, (void**)&ubs, 0);

   nlp->nunflushedvaradd = 0;

   return SCIP_OKAY;
}

/** adds variables to NLPI problem that have been added to NLP before
 * assumes that there are no unflushed variable additions or deletions (nlpFlushVarDeletions and nlpFlushVarAdditions should be called first)
 */
static
SCIP_RETCODE nlpFlushObjective(
   SCIP_NLP*             nlp,                /**< NLP data */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set                 /**< global SCIP settings */
)
{
   assert(nlp    != NULL);
   assert(blkmem != NULL);
   assert(set    != NULL);
   assert(nlp->solver != NULL);
   assert(nlp->problem != NULL);
   assert(nlp->nunflushedvaradd == 0);
   assert(nlp->nunflushedvardel == 0);
   assert(!nlp->indiving);

   if( nlp->objflushed )
      return SCIP_OKAY;

   if( nlp->objective == NULL )
   {
      /* setup SCIP objective (which is linear) */
      int*       linindices;
      SCIP_Real* lincoefs;
      SCIP_Real  coef;
      int        i, nz;

      /* assemble coefficients */
      SCIP_CALL( SCIPbufferAllocMem(set->buffer, set, (void**)&linindices, nlp->nvars_solver * sizeof(int)) );
      SCIP_CALL( SCIPbufferAllocMem(set->buffer, set, (void**)&lincoefs,   nlp->nvars_solver * sizeof(SCIP_Real)) );

      nz = 0;
      for( i = 0; i < nlp->nvars_solver; ++i )
      {
         assert(nlp->varmap_nlpi2nlp[i] >= 0); /* there should be no variable deletions pending */

         coef = SCIPvarGetObj(nlp->vars[nlp->varmap_nlpi2nlp[i]]);
         if( SCIPsetIsZero(set, coef) )
            continue;

         linindices[nz] = i;
         lincoefs[nz]   = coef;
         ++nz;
      }

      SCIP_CALL( SCIPnlpiSetObjective(nlp->solver, nlp->problem,
         nz, linindices, lincoefs,
         0, NULL,
         NULL, NULL,
         0.0) ); /* @todo would be nice to put SCIPgetTransObjOffset(scip) here */

      SCIPbufferFreeMem(set->buffer, (void**)&linindices, 0);
      SCIPbufferFreeMem(set->buffer, (void**)&lincoefs,   0);
   }
   else
   {
      SCIP_QUADELEM* quadelems;

      /* set user given objective */
      SCIP_CALL( nlpSetupNlpiIndices(nlp, blkmem, set, nlp->objective) );
      assert(nlp->objective->linvarsnlpiidx  == NULL || nlp->objective->nlinvars  == 0);
      assert(nlp->objective->quadvarsnlpiidx == NULL || nlp->objective->nquadvars == 0);
      assert(nlp->objective->exprtreenlpiidx == NULL || nlp->objective->exprtree  == NULL);

      if( nlp->objective->nquadelems > 0 )
      {
         int j;
         SCIP_CALL( SCIPbufferAllocMem(set->buffer, set, (void**)&quadelems, nlp->objective->nquadelems * sizeof(SCIP_QUADELEM)) );
         for( j = 0; j < nlp->objective->nquadelems; ++j )
         {
            quadelems[j].idx1 = MIN(nlp->objective->quadvarsnlpiidx[nlp->objective->quadelems[j].idx1], nlp->objective->quadvarsnlpiidx[nlp->objective->quadelems[j].idx2]);
            quadelems[j].idx2 = MAX(nlp->objective->quadvarsnlpiidx[nlp->objective->quadelems[j].idx1], nlp->objective->quadvarsnlpiidx[nlp->objective->quadelems[j].idx2]);
            quadelems[j].coef = nlp->objective->quadelems[j].coef;
         }
      }
      else
         quadelems = NULL;

      SCIP_CALL( SCIPnlpiSetObjective(nlp->solver, nlp->problem,
         nlp->objective->nlinvars, nlp->objective->linvarsnlpiidx, nlp->objective->lincoefs,
         nlp->objective->nquadelems, quadelems,
         nlp->objective->exprtreenlpiidx, nlp->objective->exprtree,
         0.0) );

      if( quadelems != NULL )
         SCIPbufferFreeMem(set->buffer, (void**)&quadelems, 0);
   }

   nlp->objflushed = TRUE;

   return SCIP_OKAY;
}

/** solves the NLP, assuming it has been flushed already
 * is used also to solve diving NLP */
static
SCIP_RETCODE nlpSolve(
   SCIP_NLP*             nlp,                /**< NLP data */
   BMS_BLKMEM*           blkmem,             /**< block memory buffers */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat                /**< problem statistics */
   )
{
   int i;

   assert(nlp    != NULL);
   assert(blkmem != NULL);
   assert(set    != NULL);
   assert(stat   != NULL);

   /* set initial guess, if available */
   if( nlp->haveinitguess )
   { /* @todo should we not set it if we had set it already? (initguessflushed...) */
      SCIP_Real* initialguess_solver;
      int nlpidx;

      assert(nlp->initialguess != NULL);

      SCIP_CALL( SCIPbufferAllocMem(set->buffer, set, (void**)&initialguess_solver, nlp->nvars_solver * sizeof(SCIP_Real)) );

      for( i = 0; i < nlp->nvars_solver; ++i )
      {
         nlpidx = nlp->varmap_nlpi2nlp[i];
         assert(nlpidx >= 0);
         assert(nlpidx < nlp->nvars);

         initialguess_solver[i] = nlp->initialguess[nlpidx];
      }
      SCIP_CALL( SCIPnlpiSetInitialGuess(nlp->solver, nlp->problem, initialguess_solver) );

      SCIPbufferFreeMem(set->buffer, (void**)&initialguess_solver, 0);
   }

   /* let NLP solver do his work */
   SCIPclockStart(stat->nlpsoltime, set);

   SCIP_CALL( SCIPnlpiSolve(nlp->solver, nlp->problem) );

   SCIPclockStop(stat->nlpsoltime, set);
   ++stat->nnlps;

   nlp->termstat = SCIPnlpiGetTermstat(nlp->solver, nlp->problem);
   nlp->solstat  = SCIPnlpiGetSolstat(nlp->solver, nlp->problem);
   switch( nlp->solstat )
   {
      case SCIP_NLPSOLSTAT_GLOBOPT:
      case SCIP_NLPSOLSTAT_LOCOPT:
      case SCIP_NLPSOLSTAT_FEASIBLE:
      case SCIP_NLPSOLSTAT_LOCINFEASIBLE:
      {
         SCIP_Real* solversol;

         /* store solution */
         SCIP_CALL( SCIPnlpiGetSolution(nlp->solver, nlp->problem, &solversol) );
         assert(solversol != NULL);

         if( nlp->primalsolution == NULL )
         {
            SCIP_ALLOC( BMSallocBlockMemoryArray(blkmem, &nlp->primalsolution, nlp->nvars) );
         }

         /* evaluate objective function */
         if( nlp->indiving && nlp->divingobj != NULL )
         {
            /* evaluate modified diving objective */
            SCIP_CALL( SCIPnlrowGetNLPActivity(nlp->divingobj, set, stat, nlp, &nlp->primalsolobjval) );
         }
         else if( nlp->objective == NULL )
         {
            /* evaluate default SCIP objective */
            nlp->primalsolobjval = 0.0;
            for( i = 0; i < nlp->nvars; ++i )
            {
               nlp->primalsolution[i] = solversol[nlp->varmap_nlp2nlpi[i]];
               nlp->primalsolobjval += SCIPvarGetObj(nlp->vars[i]) * nlp->primalsolution[i];
            }
         }
         else
         {
            /* evaluate non-default objective function */
            SCIP_CALL( SCIPnlrowGetNLPActivity(nlp->objective, set, stat, nlp, &nlp->primalsolobjval) );
         }
         break;
      }
      default:
         nlp->primalsolobjval = SCIP_INVALID;
         break;
   }

   return SCIP_OKAY;
}

/** event handling for variable events */
static
SCIP_DECL_EVENTEXEC(eventExecNlp)
{
   SCIP_EVENTTYPE etype;
   SCIP_VAR*      var;

   assert(scip      != NULL);
   assert(eventhdlr != NULL);
   assert(event     != NULL);
   assert(eventdata != NULL);

   assert((SCIP_NLP*)eventdata == scip->nlp);

   etype = SCIPeventGetType(event);
   var   = SCIPeventGetVar(event);

   if( SCIP_EVENTTYPE_VARADDED & etype )
   {
      SCIPdebugMessage( "-> handling varadd event, variable <%s>\n", SCIPvarGetName(var) );
      SCIP_CALL( SCIPnlpAddVar(scip->nlp, SCIPblkmem(scip), scip->set, var) );
   }
   else if( SCIP_EVENTTYPE_VARDELETED & etype )
   {
      SCIPdebugMessage( "-> handling vardel event, variable <%s>\n", SCIPvarGetName(var) );
      SCIP_CALL( SCIPnlpDelVar(scip->nlp, SCIPblkmem(scip), scip->set, scip->lp, var) );
   }
   else if( SCIP_EVENTTYPE_VARFIXED & etype )
   {
      SCIPdebugMessage( "-> handling variable fixation event, variable <%s>\n", SCIPvarGetName(var) );
      SCIPerrorMessage("handling of var fixation event not implemented yet\n");
      return SCIP_ERROR;
      /* variable was fixed, aggregated, or multiaggregated */
      /* @todo need to propagate this into rows!! */
   }
   else if( SCIP_EVENTTYPE_BOUNDCHANGED & etype )
   {
      SCIPdebugMessage( "-> handling bound changed event %d, variable <%s>\n", etype, SCIPvarGetName(var) );
      SCIP_CALL( nlpUpdateVarBounds(scip->nlp, var) );
   }
   else if( SCIP_EVENTTYPE_OBJCHANGED & etype )
   {
      SCIPdebugMessage( "-> handling objchg event, variable <%s>\n", SCIPvarGetName(var) );
      SCIP_CALL( nlpUpdateScipObjCoef(scip->nlp, var) );
   }
   else
   {
      SCIPerrorMessage("unexpected event %d on variable <%s>\n", etype, SCIPvarGetName(var) );
      return SCIP_ERROR;
   }

   return SCIP_OKAY;
}

/*
 * public NLP methods
 */

/** includes NLP specific plugins (e.g., event handler) and parameters */
SCIP_RETCODE SCIPnlpInclude(
   SCIP_SET*             set,                /**< global SCIP settings */
   BMS_BLKMEM*           blkmem              /**< block memory */
)
{
   SCIP_EVENTHDLR* eventhdlr;

   assert(set != NULL);
   assert(set->stage == SCIP_STAGE_INIT);

   /* check whether event handler is already present */
   if( SCIPsetFindEventhdlr(set, EVENTHDLR_NAME) != NULL )
   {
      SCIPerrorMessage("event handler <" EVENTHDLR_NAME "> already included.\n");
      return SCIP_INVALIDDATA;
   }

   SCIP_CALL( SCIPeventhdlrCreate(&eventhdlr, EVENTHDLR_NAME, EVENTHDLR_DESC,
         NULL, NULL, NULL, NULL, NULL, NULL, NULL, eventExecNlp, NULL) );
   SCIP_CALL( SCIPsetIncludeEventhdlr(set, eventhdlr) );

   return SCIP_OKAY;
}

/** construct a new empty NLP */
SCIP_RETCODE SCIPnlpCreate(
   SCIP_NLP**            nlp,                /**< NLP handler, call by reference */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   const char*           name,               /**< problem name */
   int                   nvars_estimate      /**< an estimate on the number of variables that may be added to the NLP later */
   )
{
   assert(nlp  != NULL);
   assert(blkmem != NULL);
   assert(set  != NULL);
   assert(stat != NULL);
   assert(name != NULL);

   SCIP_ALLOC( BMSallocMemory(nlp) );

   /* NLP solver and problem */
   assert(set->nlp_solver != NULL);
   if( set->nlp_solver[0] == '\0' )
   { /* take solver with highest priority */
      if( set->nnlpis == 0 )
      {
         SCIPerrorMessage("Attempt to create an NLP, but no NLP solver plugin (NLPI) is available.\n");
         return SCIP_ERROR;
      }
      else
      {
         assert(set->nlpis != NULL);
         (*nlp)->solver = set->nlpis[set->nnlpis-1];
      }
   }
   else
   { /* find user specified NLP solver */
      (*nlp)->solver = SCIPsetFindNlpi(set, set->nlp_solver);
      if( (*nlp)->solver == NULL )
      {
         SCIPerrorMessage("Selected NLP solver <%s> not available.\n", set->nlp_solver);
         return SCIP_PLUGINNOTFOUND;
      }
   }
   assert((*nlp)->solver != NULL);
   SCIP_CALL( SCIPnlpiCreateProblem((*nlp)->solver, &(*nlp)->problem, "scip_nlp") );

   /* status */
   (*nlp)->nunflushedvaradd   = 0;
   (*nlp)->nunflushedvardel   = 0;
   (*nlp)->nunflushednlrowadd = 0;
   (*nlp)->nunflushednlrowdel = 0;
   (*nlp)->isrelax    = TRUE;
   (*nlp)->isconvex   = TRUE;
   (*nlp)->indiving   = FALSE;

   /* variables in problem and NLPI problem */
   (*nlp)->nvars = 0;
   (*nlp)->sizevars = 0;
   (*nlp)->vars = NULL;
   SCIP_CALL( SCIPhashmapCreate(&(*nlp)->varhash, blkmem, SCIPcalcHashtableSize(nvars_estimate)) );

   (*nlp)->nvars_solver = 0;
   (*nlp)->sizevars_solver = 0;
   (*nlp)->varmap_nlp2nlpi = NULL;
   (*nlp)->varmap_nlpi2nlp = NULL;

   /* nonlinear rows in problem and NLPI problem */
   (*nlp)->nnlrows = 0;
   (*nlp)->sizenlrows = 0;
   (*nlp)->nlrows = NULL;

   (*nlp)->nnlrows_solver = 0;
   (*nlp)->sizenlrows_solver = 0;
   (*nlp)->nlrowmap_nlpi2nlp = NULL;

   /* objective function */
   (*nlp)->objective = NULL;
   (*nlp)->objflushed = TRUE;
   (*nlp)->divingobj = NULL;

   /* initial guess */
   (*nlp)->haveinitguess = FALSE;
   (*nlp)->initialguess = NULL;

   /* solution of NLP */
   (*nlp)->primalsolution  = NULL;
   (*nlp)->primalsolobjval = SCIP_INVALID;
   (*nlp)->solstat         = SCIP_NLPSOLSTAT_UNKNOWN;
   (*nlp)->termstat        = SCIP_NLPTERMSTAT_OTHER;

   /* event handling: catch variable addition and deletion events */
   (*nlp)->eventhdlr = SCIPsetFindEventhdlr(set, EVENTHDLR_NAME);
   if( (*nlp)->eventhdlr == NULL )
   {
      SCIPerrorMessage("NLP eventhandler <"EVENTHDLR_NAME"> not found.\n");
      return SCIP_PLUGINNOTFOUND;
   }
   SCIP_CALL( SCIPeventfilterAdd(set->scip->eventfilter, blkmem, set,
      SCIP_EVENTTYPE_VARADDED | SCIP_EVENTTYPE_VARDELETED,
      (*nlp)->eventhdlr, (SCIP_EVENTDATA*)(*nlp), &(*nlp)->globalfilterpos) );

   /* miscellaneous */
   SCIP_ALLOC( BMSduplicateBlockMemoryArray(blkmem, &(*nlp)->name, name, strlen(name)+1) );

   return SCIP_OKAY;
}

/** frees NLP data object */
SCIP_RETCODE SCIPnlpFree(
   SCIP_NLP**            nlp,                /**< pointer to NLP data object */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_LP*              lp                  /**< SCIP LP, needed for releasing variables */
   )
{
   assert(nlp    != NULL);
   assert(*nlp   != NULL);
   assert(blkmem != NULL);
   assert(set    != NULL);

   /* drop global events (variable addition and deletion) */
   SCIP_CALL( SCIPeventfilterDel(set->scip->eventfilter, blkmem, set,
      SCIP_EVENTTYPE_VARADDED | SCIP_EVENTTYPE_VARDELETED,
      (*nlp)->eventhdlr, (SCIP_EVENTDATA*)(*nlp), (*nlp)->globalfilterpos) );

   SCIP_CALL( SCIPnlpReset(*nlp, blkmem, set, lp) );
   assert((*nlp)->objective == NULL);
   assert((*nlp)->nnlrows == 0);
   assert((*nlp)->nnlrows_solver == 0);
   assert((*nlp)->nvars == 0);
   assert((*nlp)->nvars_solver == 0);
   assert((*nlp)->primalsolution == NULL);
   assert((*nlp)->initialguess == NULL);

   BMSfreeBlockMemoryArray(blkmem, &(*nlp)->name, strlen((*nlp)->name)+1);

   /* free nonlinear rows arrays */
   BMSfreeBlockMemoryArrayNull(blkmem, &(*nlp)->nlrowmap_nlpi2nlp, (*nlp)->sizenlrows_solver);
   BMSfreeBlockMemoryArrayNull(blkmem, &(*nlp)->nlrows, (*nlp)->sizenlrows);

   /* free variables arrays */
   BMSfreeBlockMemoryArrayNull(blkmem, &(*nlp)->varmap_nlp2nlpi, (*nlp)->sizevars);
   BMSfreeBlockMemoryArrayNull(blkmem, &(*nlp)->varmap_nlpi2nlp, (*nlp)->sizevars_solver);
   SCIPhashmapFree(&(*nlp)->varhash);
   BMSfreeBlockMemoryArrayNull(blkmem, &(*nlp)->vars, (*nlp)->sizevars);

   /* free NLPI problem */
   SCIP_CALL( SCIPnlpiFreeProblem((*nlp)->solver, &(*nlp)->problem) );

   /* free NLP data structure */
   BMSfreeMemory(nlp);

   return SCIP_OKAY;
}

/** resets the NLP to the empty NLP by removing all variables and rows from NLP,
 *  releasing all rows, and flushing the changes to the NLP solver
 */
SCIP_RETCODE SCIPnlpReset(
   SCIP_NLP*             nlp,                /**< NLP data */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_LP*              lp                  /**< SCIP LP, needed for releasing variables */
   )
{
   int i;

   assert(nlp    != NULL);
   assert(blkmem != NULL);
   assert(set    != NULL);

   if( nlp->indiving )
   {
      SCIP_CALL( SCIPnlpEndDive(nlp, blkmem, set) );
   }

   BMSfreeBlockMemoryArrayNull(blkmem, &nlp->initialguess, nlp->nvars);
   nlp->haveinitguess = FALSE;

   BMSfreeBlockMemoryArrayNull(blkmem, &nlp->primalsolution, nlp->nvars);

   SCIP_CALL( SCIPnlpSetObjective(nlp, blkmem, set, NULL) );

   for(i = nlp->nnlrows - 1; i >= 0; --i)
   {
      SCIP_CALL( nlpDelNlRowPos(nlp, blkmem, set, i) );
   }

   for(i = nlp->nvars - 1; i >= 0; --i)
   {
      SCIP_CALL( nlpDelVarPos(nlp, blkmem, set, lp, i) );
   }

   SCIP_CALL( SCIPnlpFlush(nlp, blkmem, set) );

   return SCIP_OKAY;
}

/** currently a dummy function that always returns TRUE */
SCIP_Bool SCIPnlpHasCurrentNodeNLP(
   SCIP_NLP*             nlp                 /**< NLP data */
)
{
   return TRUE;
}

/** ensures, that variables array of NLP can store at least num entries */
SCIP_RETCODE SCIPnlpEnsureVarsSize(
   SCIP_NLP*             nlp,                /**< NLP data */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   int                   num                 /**< minimum number of entries to store */
   )
{
   assert(nlp    != NULL);
   assert(blkmem != NULL);
   assert(set    != NULL);
   assert(nlp->nvars <= nlp->sizevars);

   if( num > nlp->sizevars )
   {
      int newsize;

      newsize = SCIPsetCalcMemGrowSize(set, num);
      SCIP_ALLOC( BMSreallocBlockMemoryArray(blkmem, &nlp->vars,            nlp->sizevars, newsize) );
      SCIP_ALLOC( BMSreallocBlockMemoryArray(blkmem, &nlp->varmap_nlp2nlpi, nlp->sizevars, newsize) );
      if( nlp->initialguess != NULL )
      {
         SCIP_ALLOC( BMSreallocBlockMemoryArray(blkmem, &nlp->initialguess, nlp->sizevars, newsize) );
      }
      if( nlp->primalsolution != NULL )
      {
         SCIP_ALLOC( BMSreallocBlockMemoryArray(blkmem, &nlp->primalsolution, nlp->sizevars, newsize) );
      }

      nlp->sizevars = newsize;
   }
   assert(num <= nlp->sizevars);

   return SCIP_OKAY;
}

/** adds a variable to the NLP and captures the variable */
SCIP_RETCODE SCIPnlpAddVar(
   SCIP_NLP*             nlp,                /**< NLP data */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_VAR*             var                 /**< variable */
   )
{
   assert(nlp != NULL);
   assert(blkmem != NULL);
   assert(set != NULL);
   assert(var != NULL);
   assert(SCIPvarIsTransformed(var));
   assert(!SCIPhashmapExists(nlp->varhash, var));

   if( nlp->indiving )
   {
      SCIPerrorMessage("cannot add variable during NLP diving\n");
      return SCIP_ERROR;
   }

   SCIP_CALL( nlpAddVars(nlp, blkmem, set, 1, &var) );

   return SCIP_OKAY;
}

/** adds a set of variables to the NLP and captures the variables */
SCIP_RETCODE SCIPnlpAddVars(
   SCIP_NLP*             nlp,                /**< NLP data */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   int                   nvars,              /**< number of variables to add */
   SCIP_VAR**            vars                /**< variables to add */
   )
{
   assert(nlp != NULL);
   assert(blkmem != NULL);
   assert(set != NULL);
   assert(vars != NULL || nvars == 0);

   if( nlp->indiving && nvars > 0)
   {
      SCIPerrorMessage("cannot add variables during NLP diving\n");
      return SCIP_ERROR;
   }

   SCIP_CALL( nlpAddVars(nlp, blkmem, set, nvars, vars) );

   return SCIP_OKAY;
}

/** deletes a variable from the NLP and releases the variable */
SCIP_RETCODE SCIPnlpDelVar(
   SCIP_NLP*             nlp,                /**< NLP data */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_LP*              lp,                 /**< SCIP LP, needed to release variable */
   SCIP_VAR*             var                 /**< variable */
   )
{
   int varpos;

   assert(nlp    != NULL);
   assert(blkmem != NULL);
   assert(set    != NULL);
   assert(var    != NULL);

   if( !SCIPhashmapExists(nlp->varhash, var) )
   {
      SCIPerrorMessage("variable <%s> not found in NLP, cannot delete\n", SCIPvarGetName(var));
      return SCIP_ERROR;
   }

   if( nlp->indiving )
   {
      SCIPerrorMessage("cannot delete variable during NLP diving\n");
      return SCIP_ERROR;
   }

   varpos = (int) (size_t) SCIPhashmapGetImage(nlp->varhash, var);

   SCIP_CALL( nlpDelVarPos(nlp, blkmem, set, lp, varpos) );

   return SCIP_OKAY;
}


/** ensures, that nonlinear rows array of NLP can store at least num entries */
SCIP_RETCODE SCIPnlpEnsureNlRowsSize(
   SCIP_NLP*             nlp,                /**< NLP data */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   int                   num                 /**< minimum number of entries to store */
   )
{
   assert(nlp    != NULL);
   assert(blkmem != NULL);
   assert(set    != NULL);
   assert(nlp->nnlrows <= nlp->sizenlrows);

   if( num > nlp->sizenlrows )
   {
      int newsize;

      newsize = SCIPsetCalcMemGrowSize(set, num);
      SCIP_ALLOC( BMSreallocBlockMemoryArray(blkmem, &nlp->nlrows, nlp->sizenlrows, newsize) );

      nlp->sizenlrows = newsize;
   }
   assert(num <= nlp->sizenlrows);

   return SCIP_OKAY;
}

/** adds a nonlinear row to the NLP and captures it
 * all variables of the row need to be present in the NLP */
SCIP_RETCODE SCIPnlpAddNlRow(
   SCIP_NLP*             nlp,                /**< NLP data */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics data */
   SCIP_NLROW*           nlrow               /**< nonlinear row */
   )
{
   assert(nlp   != NULL);
   assert(nlrow != NULL);

   if( nlp->indiving )
   {
      SCIPerrorMessage("cannot add row during NLP diving\n");
      return SCIP_ERROR;
   }

   SCIP_CALL( nlpAddNlRows(nlp, blkmem, set, stat, 1, &nlrow) );

   return SCIP_OKAY;
}

/** adds nonlinear rows to the NLP and captures them
 * all variables of the row need to be present in the NLP */
SCIP_RETCODE SCIPnlpAddNlRows(
   SCIP_NLP*             nlp,                /**< NLP data */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics data */
   int                   nnlrows,            /**< number of rows to add */
   SCIP_NLROW**          nlrows              /**< rows to add */
   )
{
   assert(nlp    != NULL);
   assert(nlrows != NULL || nnlrows == 0);

   if( nnlrows == 0 )
      return SCIP_OKAY;

   if( nlp->indiving )
   {
      SCIPerrorMessage("cannot add rows during NLP diving\n");
      return SCIP_ERROR;
   }

   SCIP_CALL( nlpAddNlRows(nlp, blkmem, set, stat, nnlrows, nlrows) );

   return SCIP_OKAY;
}

/** deletes a nonlinear row from the NLP
 * does nothing if nonlinear row is not in NLP */
SCIP_RETCODE SCIPnlpDelNlRow(
   SCIP_NLP*             nlp,                /**< NLP data */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_NLROW*           nlrow               /**< nonlinear row */
   )
{
   assert(nlp    != NULL);
   assert(blkmem != NULL);
   assert(set    != NULL);
   assert(nlrow  != NULL);

   /* if row not in NLP, nothing to do */
   if( nlrow->nlpindex == -2 )
      return SCIP_OKAY;

   if( nlrow->nlpindex == -1 )
   {
      SCIPerrorMessage("cannot remove objective function by using SCIPnlpDelNlRow\n");
      return SCIP_ERROR;
   }
   assert(nlrow->nlpindex >= 0);
   assert(nlrow->nlpindex < nlp->nnlrows);

   if( nlp->indiving )
   {
      SCIPerrorMessage("cannot delete row during NLP diving\n");
      return SCIP_ERROR;
   }

   SCIP_CALL( nlpDelNlRowPos(nlp, blkmem, set, nlrow->nlpindex) );

   return SCIP_OKAY;
}

/** sets the objective function
 * If a nonliner row is given, then the row function is used as objective function and its bounds are ignored.
 * The row is captured.
 * If NULL is given, then a linear objective with coefficients taken from the SCIP problem is used (i.e., objective coefficients as stored in variables that are part of the NLP).
 */
SCIP_RETCODE SCIPnlpSetObjective(
   SCIP_NLP*             nlp,                /**< NLP data */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_NLROW*           obj                 /**< new objective as nonlinear row, or NULL for SCIP objective */
   )
{
#ifndef NDEBUG
   int i;
#endif

   assert(nlp    != NULL);
   assert(blkmem != NULL);
   assert(set    != NULL);

   if( nlp->indiving )
   {
      SCIPerrorMessage("cannot set objective during NLP diving (use SCIPchgVarObjDiveNLP to change single coefficients)\n");
      return SCIP_ERROR;
   }

   /* if previous and new objective are the same, nothing to do */
   if( nlp->objective == obj )
      return SCIP_OKAY;

   /* release previous objective, if present */
   if( nlp->objective != NULL )
   {
      /* this row is not in the NLP and NLPI anymore from now on */
      nlp->objective->nlpindex = -2;
      nlp->objective->nlpiindex = -2;
      SCIP_CALL( SCIPnlrowRelease(&nlp->objective, blkmem, set) );
      assert(nlp->objective == NULL);
   }

   /* install the new objective function */
   if( obj != NULL )
   {
      nlp->objective = obj;
      obj->nlpindex = -1;
      SCIPnlrowCapture(obj);

#ifndef NDEBUG
      /* assert that variables of row are in NLP */
      for( i = 0; i < obj->nlinvars; ++i )
         assert(SCIPhashmapExists(nlp->varhash, obj->linvars[i]));

      for( i = 0; i < obj->nquadvars; ++i )
         assert(SCIPhashmapExists(nlp->varhash, obj->quadvars[i]));

      if( obj->exprtree )
      {
         int n;

         n = SCIPexprtreeGetNVars(obj->exprtree);
         assert(SCIPexprtreeGetVars(obj->exprtree) != NULL || n == 0);

         for( i = 0; i < n; ++i )
            assert(SCIPhashmapExists(nlp->varhash, SCIPexprtreeGetVars(obj->exprtree)[i]));
      }
#endif
   }

   nlp->objflushed = FALSE;

   /* if we were feasible before, then we stay feasible
    * if we were locally or globally optimal, then we are now still feasible
    * if we were infeasible, then we are still infeasible
    * if we were unbounded, then we may not be unbounded anymore
    */
   if( nlp->solstat <= SCIP_NLPSOLSTAT_LOCOPT )
      nlp->solstat = SCIP_NLPSOLSTAT_FEASIBLE;
   else if( nlp->solstat == SCIP_NLPSOLSTAT_UNBOUNDED )
      nlp->solstat = SCIP_NLPSOLSTAT_UNKNOWN;

   return SCIP_OKAY;
}

/** applies all cached changes to the NLP solver */
SCIP_RETCODE SCIPnlpFlush(
   SCIP_NLP*             nlp,                /**< current NLP data */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set                 /**< global SCIP settings */
   )
{
   assert(nlp    != NULL);
   assert(blkmem != NULL);
   assert(set    != NULL);

   if( nlp->indiving )
   {
      SCIPerrorMessage("cannot flush NLP during NLP diving\n");
      return SCIP_ERROR;
   }

   /* flush removals of nonlinear rows and variables */
   SCIP_CALL( nlpFlushNlRowDeletions(nlp, blkmem, set) );
   SCIP_CALL( nlpFlushVarDeletions(nlp, blkmem, set) );
   assert(nlp->nunflushednlrowdel == 0);
   assert(nlp->nunflushedvardel   == 0);

   /* flush addition of variables, setting of objective, and addition of rows */
   SCIP_CALL( nlpFlushVarAdditions(nlp, blkmem, set) );
   SCIP_CALL( nlpFlushObjective(nlp, blkmem, set) );
   SCIP_CALL( nlpFlushNlRowAdditions(nlp, blkmem, set) );
   assert(nlp->nunflushedvaradd == 0);
   assert(nlp->objflushed == TRUE);
   assert(nlp->nunflushednlrowadd == 0);

   assert(nlp->nvars   == nlp->nvars_solver);
   assert(nlp->nnlrows == nlp->nnlrows_solver);

   return SCIP_OKAY;
}

/** solves the NLP */
SCIP_RETCODE SCIPnlpSolve(
   SCIP_NLP*             nlp,                /**< NLP data */
   BMS_BLKMEM*           blkmem,             /**< block memory buffers */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat                /**< problem statistics */
   )
{
   assert(nlp    != NULL);
   assert(blkmem != NULL);
   assert(set    != NULL);
   assert(stat   != NULL);

   if( nlp->indiving )
   {
      SCIPerrorMessage("cannot solve NLP during NLP diving (use SCIPsolveDiveNLP)\n");
      return SCIP_ERROR;
   }

   SCIP_CALL( SCIPnlpFlush(nlp, blkmem, set) );

   SCIP_CALL( nlpSolve(nlp, blkmem, set, stat) );

   return SCIP_OKAY;
}

/** gets objective value of current NLP */
SCIP_Real SCIPnlpGetObjval(
   SCIP_NLP*             nlp                 /**< current NLP data */
   )
{
   assert(nlp != NULL);

   return nlp->primalsolobjval;
}

/** gives current pseudo objective value */
SCIP_RETCODE SCIPnlpGetPseudoObjval(
   SCIP_NLP*             nlp,                /**< current NLP data */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_Real*            pseudoobjval        /**< buffer to store pseudo objective value */
   )
{
   assert(nlp != NULL);
   assert(pseudoobjval != NULL);

   if( nlp->divingobj != NULL )
   {
      assert(nlp->indiving);
      SCIP_CALL( SCIPnlrowGetPseudoActivity(nlp->divingobj, set, stat, pseudoobjval) );
   }
   else if( nlp->objective == NULL )
   {
      int i;

      *pseudoobjval = 0.0;
      for( i = 0; i < nlp->nvars; ++i )
         *pseudoobjval += SCIPvarGetObj(nlp->vars[i]) * SCIPvarGetBestBound(nlp->vars[i]);
   }
   else
   {
      SCIP_CALL( SCIPnlrowGetPseudoActivity(nlp->objective, set, stat, pseudoobjval) );
   }

   return SCIP_OKAY;
}

/** provides current primal solution in new SCIP_SOL data structure
 * *sol is set to NULL if no NLP solution is available */
SCIP_RETCODE SCIPnlpGetSol(
   SCIP_NLP*             nlp,                /**< NLP data structure */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_PRIMAL*          primal,             /**< primal data */
   SCIP_TREE*            tree,               /**< branch and bound tree */
   SCIP_SOL**            sol,                /**< buffer where to store pointer to new solution */
   SCIP_HEUR*            heur                /**< heuristic that solved the NLP, or NULL if not from a heuristic */
)
{
   int i;

   assert(nlp != NULL);
   assert(sol != NULL);
   assert(nlp->solstat > SCIP_NLPSOLSTAT_FEASIBLE || nlp->primalsolution != NULL);

   if( nlp->solstat > SCIP_NLPSOLSTAT_FEASIBLE || nlp->primalsolution == NULL )
   {
      *sol = NULL;
      return SCIP_OKAY;
   }

   SCIP_CALL( SCIPsolCreate(sol, blkmem, set, stat, primal, tree, heur) );
   assert(*sol != NULL);

   for( i = 0; i < nlp->nvars; ++i )
   {
      SCIP_CALL( SCIPsolSetVal(*sol, set, stat, tree, nlp->vars[i], nlp->primalsolution[i]) );
   }

   return SCIP_OKAY;
}

/** removes all redundant nonlinear rows */
SCIP_RETCODE SCIPnlpRemoveRedundantNlRows(
   SCIP_NLP*             nlp,                /**< current NLP data */
   BMS_BLKMEM*           blkmem,             /**< block memory buffers */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat                /**< problem statistics */
   )
{
   SCIP_NLPSOLSTAT solstatus;
   SCIP_Bool isredundant;
   int i;

   assert(nlp    != NULL);
   assert(blkmem != NULL);
   assert(set    != NULL);
   assert(stat   != NULL);

   if( nlp->nnlrows == 0 )
      return SCIP_OKAY;

   if( nlp->indiving )
   {
      SCIPerrorMessage("cannot remove redundant rows during NLP diving\n");
      return SCIP_ERROR;
   }

   /* removing redundant rows should not change the solution status, so we reset it at the end */
   solstatus = nlp->solstat;

   for( i = 0; i < nlp->nnlrows; ++i )
   {
      SCIP_CALL( SCIPnlrowIsRedundant(nlp->nlrows[i], set, stat, &isredundant) );
      if( isredundant )
      {
         SCIP_CALL( nlpDelNlRowPos(nlp, blkmem, set, i) );
      }
   }

   return SCIP_OKAY;
}

/** set initial guess (approximate primal solution) for next solve
 * array initguess must be NULL or have length at least SCIPnlpGetNVars */
SCIP_RETCODE SCIPnlpSetInitialGuess(
   SCIP_NLP*             nlp,                /**< current NLP data */
   BMS_BLKMEM*           blkmem,             /**< block memory buffers */
   SCIP_Real*            initguess           /**< new initial guess, or NULL to clear previous one */
   )
{
   assert(nlp    != NULL);
   assert(blkmem != NULL);
   assert(nlp->solver  != NULL);
   assert(nlp->problem != NULL);

   /* if user wants to let NLP solver choose start point, then invalidate current initial guess both in NLP and in NLPI */
   if( initguess == NULL )
   {
      nlp->haveinitguess = FALSE;
      SCIP_CALL( SCIPnlpiSetInitialGuess(nlp->solver, nlp->problem, NULL) );
   }

   if( nlp->initialguess != NULL )
   {
      BMScopyMemoryArray(nlp->initialguess, initguess, nlp->nvars);
   }
   else
   {
      SCIP_ALLOC( BMSduplicateBlockMemoryArray(blkmem, &nlp->initialguess, initguess, nlp->nvars) );
   }
   nlp->haveinitguess = TRUE;

   return SCIP_OKAY;
}

/** writes NLP to a file */
SCIP_RETCODE SCIPnlpWrite(
   SCIP_NLP*             nlp,                /**< current NLP data */
   SCIP_SET*             set,                /**< global SCIP settings */
   const char*           fname               /**< file name */
   )
{
   FILE* file;
   int i;

   assert(nlp != NULL);

   if( fname != NULL )
   {
      file = fopen(fname, "w");
      if( file == NULL )
      {
         SCIPerrorMessage("could not open file <%s> for writing\n", fname);
         return SCIP_ERROR;
      }
   }
   else
      file = stdout;

   SCIPmessageFPrintInfo(file, "STATISTICS\n");
   SCIPmessageFPrintInfo(file, "  NLP name: %s\n", nlp->name);
   SCIPmessageFPrintInfo(file, "  Variables: %d\n", nlp->nvars);
   SCIPmessageFPrintInfo(file, "  Rows: %d\n", nlp->nnlrows);

   SCIPmessageFPrintInfo(file, "VARIABLES\n");
   for( i = 0; i < nlp->nvars; ++i )
   {
      SCIPvarPrint(nlp->vars[i], set, file);
   }

   if( nlp->objective != NULL )
   {
      SCIPmessageFPrintInfo(file, "OBJECTIVE\n");
      SCIPnlrowPrint(nlp->objective, file);
   }

   SCIPmessageFPrintInfo(file, "NONLINEAR ROWS\n");
   for( i = 0; i < nlp->nnlrows; ++i )
   {
      SCIPmessageFPrintInfo(file, "  ");
      SCIPnlrowPrint(nlp->nlrows[i], file);
   }

   if( fname != NULL )
   {
      fclose(file);
   }

   return SCIP_OKAY;
}

/*
 * NLP diving methods
 */

/** signals start of diving */
SCIP_RETCODE SCIPnlpStartDive(
   SCIP_NLP*             nlp,                /**< current NLP data */
   BMS_BLKMEM*           blkmem,             /**< block memory buffers */
   SCIP_SET*             set                 /**< global SCIP settings */
   )
{
   assert(nlp != NULL);

   if( nlp->indiving )
   {
      SCIPerrorMessage("NLP is already in diving mode\n");
      return SCIP_ERROR;
   }

   SCIP_CALL( SCIPnlpFlush(nlp, blkmem, set) );

   nlp->indiving = TRUE;

   return SCIP_OKAY;
}

/** resets the bound and objective changes made during diving and disables diving mode */
extern
SCIP_RETCODE SCIPnlpEndDive(
   SCIP_NLP*             nlp,                /**< current NLP data */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set                 /**< global SCIP settings */
   )
{
   int i;
   int* varidx;
   SCIP_Real* varlb;
   SCIP_Real* varub;

   assert(nlp != NULL);
   assert(set != NULL);
   assert(nlp->nvars == nlp->nvars_solver);

   if( !nlp->indiving )
   {
      SCIPerrorMessage("NLP not in diving mode, cannot end dive\n");
      return SCIP_ERROR;
   }

   /* reset variable bounds in NLPI problem to their current values */
   SCIP_CALL( SCIPbufferAllocMem(set->buffer, set, (void**)&varidx, nlp->nvars * sizeof(int)) );
   SCIP_CALL( SCIPbufferAllocMem(set->buffer, set, (void**)&varlb,  nlp->nvars * sizeof(SCIP_Real)) );
   SCIP_CALL( SCIPbufferAllocMem(set->buffer, set, (void**)&varub,  nlp->nvars * sizeof(SCIP_Real)) );
   for( i = 0; i < nlp->nvars; ++i )
   {
      varidx[i] = i;
      varlb[i] = SCIPvarGetLbLocal(nlp->vars[nlp->varmap_nlpi2nlp[i]]);
      varub[i] = SCIPvarGetUbLocal(nlp->vars[nlp->varmap_nlpi2nlp[i]]);
   }

   SCIP_CALL( SCIPnlpiChgVarBounds(nlp->solver, nlp->problem, nlp->nvars, varidx, varlb, varub) );

   SCIPbufferFreeMem(set->buffer, (void**)&varidx, 0);
   SCIPbufferFreeMem(set->buffer, (void**)&varlb,  0);
   SCIPbufferFreeMem(set->buffer, (void**)&varub,  0);

   /* clear diving objective, if one was used (i.e., if SCIPnlpChgVarObjDive had been called)
    * the objective in the NLPI will be reset in the next flush */
   if( nlp->divingobj != NULL )
   {
      SCIP_CALL( SCIPnlrowRelease(&nlp->divingobj, blkmem, set) );
      assert(nlp->divingobj == NULL);
      assert(nlp->objflushed == FALSE);
   }

   /* we do not have a valid solution anymore */
   nlp->solstat  = SCIP_NLPSOLSTAT_UNKNOWN;
   nlp->termstat = SCIP_NLPTERMSTAT_OTHER;
   nlp->primalsolobjval = SCIP_INVALID;

   nlp->indiving = FALSE;

   return SCIP_OKAY;
}

/** changes coefficient of variable in diving NLP */
SCIP_RETCODE SCIPnlpChgVarObjDive(
   SCIP_NLP*             nlp,                /**< current NLP data */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_VAR*             var,                /**< variable which coefficient to change */
   SCIP_Real             coef                /**< new linear coefficient of variable in objective */
   )
{
   int pos;
   int objidx;

   assert(nlp != NULL);
   assert(var != NULL);
   assert(SCIPhashmapExists(nlp->varhash, var));
   assert(nlp->indiving);

   /* get position of variable in NLPI problem */
   pos = (int) (size_t) SCIPhashmapGetImage(nlp->varhash, var);
   pos = nlp->varmap_nlp2nlpi[pos];
   assert(pos >= 0);

   /* set coefficient in NLPI problem objective */
   objidx = -1;
   SCIP_CALL( SCIPnlpiChgLinearCoefs(nlp->solver, nlp->problem, objidx, 1, &pos, &coef) );

   /* create diving objective as copy of original objective, if not done yet */
   if( nlp->divingobj == NULL )
   {
      if( nlp->objective == NULL )
      {
         SCIP_CALL( SCIPnlrowCreateCopy(&nlp->divingobj, blkmem, set, nlp->objective) );
      }
      else
      {
         /* setup nlrow corresponding to SCIP objective function */
         SCIP_Real* coefs;
         int        i;

         SCIP_CALL( SCIPbufferAllocMem(set->buffer, set, (void**)&coefs, nlp->nvars * sizeof(SCIP_Real)) );
         for( i = 0; i < nlp->nvars; ++i )
            coefs[i] = SCIPvarGetObj(nlp->vars[i]);

         SCIP_CALL( SCIPnlrowCreate(&nlp->divingobj, blkmem, set, "divingobj",
            nlp->nvars, nlp->vars, coefs,
            0, NULL, 0, NULL,
            NULL,
            -SCIPsetInfinity(set), SCIPsetInfinity(set)) );

         SCIPbufferFreeMem(set->buffer, (void**)&coefs, 0);
      }
      assert(nlp->divingobj != NULL);
   }

   /* modify coefficient in diving objective */
   SCIP_CALL( SCIPnlrowChgLinearCoef(nlp->divingobj, blkmem, set, nlp, var, coef) );

   /* remember that we have to store objective after diving ended */
   nlp->objflushed = FALSE;

   return SCIP_OKAY;
}

/** changes bounds of variable in diving NLP */
SCIP_RETCODE SCIPnlpChgVarBoundsDive(
   SCIP_NLP*             nlp,                /**< current NLP data */
   SCIP_VAR*             var,                /**< variable which coefficient to change */
   SCIP_Real             lb,                 /**< new lower bound of variable */
   SCIP_Real             ub                  /**< new upper bound of variable */
   )
{
   int pos;

   assert(nlp != NULL);
   assert(var != NULL);
   assert(SCIPhashmapExists(nlp->varhash, var));
   assert(nlp->indiving);

   /* get position of variable in NLPI problem */
   pos = (int) (size_t) SCIPhashmapGetImage(nlp->varhash, var);
   pos = nlp->varmap_nlp2nlpi[pos];
   assert(pos >= 0);

   /* set new bounds in NLPI */
   SCIP_CALL( SCIPnlpiChgVarBounds(nlp->solver, nlp->problem, 1, &pos, &lb, &ub) );

   return SCIP_OKAY;
}

/** changes bounds of a set of variables in diving NLP */
SCIP_RETCODE SCIPnlpChgVarsBoundsDive(
   SCIP_NLP*             nlp,                /**< current NLP data */
   SCIP_SET*             set,                /**< global SCIP settings */
   int                   nvars,              /**< number of variables which bounds to change */
   SCIP_VAR**            vars,               /**< variables which bounds to change */
   SCIP_Real*            lbs,                /**< new lower bounds of variables */
   SCIP_Real*            ubs                 /**< new upper bounds of variables */
   )
{
   int i;
   int* poss;

   assert(nlp  != NULL);
   assert(vars != NULL || nvars == 0);
   assert(nlp->indiving);
   assert(lbs  != NULL || nvars == 0);
   assert(ubs  != NULL || nvars == 0);

   if( nvars == 0 )
      return SCIP_OKAY;

   SCIP_CALL( SCIPbufferAllocMem(set->buffer, set, (void**)&poss, nvars * sizeof(int)) );

   for( i = 0; i < nvars; ++i )
   {
      assert(SCIPhashmapExists(nlp->varhash, vars[i]));

      /* get position of variable in NLPI problem */
      poss[i] = (int) (size_t) SCIPhashmapGetImage(nlp->varhash, vars[i]);
      poss[i] = nlp->varmap_nlp2nlpi[poss[i]];
      assert(poss[i] >= 0);
   }

   /* set new bounds in NLPI */
   SCIP_CALL( SCIPnlpiChgVarBounds(nlp->solver, nlp->problem, nvars, poss, lbs, ubs) );

   SCIPbufferFreeMem(set->buffer, (void**)&poss, 0);

   return SCIP_OKAY;
}


/** solves diving NLP */
SCIP_RETCODE SCIPnlpSolveDive(
   SCIP_NLP*             nlp,                /**< current NLP data */
   BMS_BLKMEM*           blkmem,             /**< block memory buffers */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat                /**< problem statistics */
   )
{
   SCIP_CALL( nlpSolve(nlp, blkmem, set, stat) );

   return SCIP_OKAY;
}


/*
 * public NLP methods
 */

/** sets whether the current NLP is a convex problem, i.e., all restrictions are defined by convex functions w.r.t. current bounds */
void SCIPnlpSetIsConvex(
   SCIP_NLP*             nlp,                /**< NLP data */
   SCIP_Bool             isconvex            /**< is the current NLP a convex problem? */
   )
{
   assert(nlp != NULL);

   nlp->isconvex = isconvex;
}

/** returns whether the current NLP is a convex problem, i.e., all restrictions are defined by convex functions w.r.t. current bounds */
SCIP_Bool SCIPnlpIsConvex(
   SCIP_NLP*             nlp                 /**< NLP data */
   )
{
   assert(nlp != NULL);

   return nlp->isconvex;
}

/** gets array with variables of the NLP */
SCIP_VAR** SCIPnlpGetVars(
   SCIP_NLP*             nlp                 /**< current NLP data */
   )
{
   assert(nlp != NULL);

   return nlp->vars;
}

/** gets current number of variables in NLP */
int SCIPnlpGetNVars(
   SCIP_NLP*             nlp                 /**< current NLP data */
   )
{
   assert(nlp != NULL);

   return nlp->nvars;
}

/** gets array with nonlinear rows of the NLP */
SCIP_NLROW** SCIPnlpGetNlRows(
   SCIP_NLP*             nlp                 /**< current NLP data */
   )
{
   assert(nlp != NULL);

   return nlp->nlrows;
}

/** gets current number of nonlinear rows in NLP */
int SCIPnlpGetNNlRows(
   SCIP_NLP*             nlp                 /**< current NLP data */
   )
{
   assert(nlp != NULL);

   return nlp->nnlrows;
}

/** gets objective of the NLP
 * gives NULL if SCIP objective is used */
SCIP_NLROW* SCIPnlpGetObjective(
   SCIP_NLP*             nlp                 /**< current NLP data */
   )
{
   assert(nlp != NULL);

   if( nlp->indiving && nlp->divingobj )
      return nlp->divingobj;

   return nlp->objective;
}

/** gets the NLP solver interface */
SCIP_NLPI* SCIPnlpGetNLPI(
   SCIP_NLP*             nlp                 /**< current NLP data */
   )
{
   assert(nlp != NULL);

   return nlp->solver;
}

/** gets the NLP problem in the solver interface */
SCIP_NLPIPROBLEM* SCIPnlpGetNLPIProblem(
   SCIP_NLP*             nlp                 /**< current NLP data */
   )
{
   assert(nlp != NULL);

   return nlp->problem;
}

/** gets solution status of current NLP */
SCIP_NLPSOLSTAT SCIPnlpGetSolstat(
   SCIP_NLP*             nlp                 /**< current NLP data */
   )
{
   assert(nlp != NULL);

   return nlp->solstat;
}

/** gets termination status of last NLP solve */
SCIP_NLPTERMSTAT SCIPnlpGetTermstat(
   SCIP_NLP*             nlp                 /**< current NLP data */
   )
{
   assert(nlp != NULL);

   return nlp->termstat;
}

/** gives statistics (number of iterations, solving time, ...) of last NLP solve */
SCIP_RETCODE SCIPnlpGetStatistics(
   SCIP_NLP*             nlp,                /**< pointer to NLP datastructure */
   SCIP_NLPSTATISTICS*   statistics          /**< pointer to store statistics */
)
{
   assert(nlp != NULL);
   assert(nlp->solver != NULL);
   assert(nlp->problem != NULL);
   assert(statistics != NULL);

   SCIP_CALL( SCIPnlpiGetStatistics(nlp->solver, nlp->problem, statistics) );

   return SCIP_OKAY;
}

/** indicates whether a feasible solution for the current NLP is available
 * thus, returns whether the solution status <= feasible  */
SCIP_Bool SCIPnlpHasSolution(
   SCIP_NLP*             nlp                 /**< current NLP data */
   )
{
   assert(nlp != NULL);

   return nlp->solstat <= SCIP_NLPSOLSTAT_FEASIBLE;
}

/** gets values of current primal NLP solution
 * returns NULL if no solution available
 * use SCIPnlpGetSolstat to get information on whether solution is optimal or just feasible
 * use SCIPnlpGetVars to get variables corresponding to solution values */
SCIP_Real* SCIPnlpGetSolVals(
   SCIP_NLP*             nlp                 /**< current NLP data */
   )
{
   assert(nlp != NULL);
   assert(nlp->solstat > SCIP_NLPSOLSTAT_FEASIBLE || nlp->primalsolution != NULL);

   if( nlp->solstat > SCIP_NLPSOLSTAT_FEASIBLE || nlp->primalsolution == NULL )
      return NULL;

   return nlp->primalsolution;
}

/** gets value of variable in current NLP solution
 * returns SCIP_INVALID if no solution available */
SCIP_Real SCIPnlpGetVarSolVal(
   SCIP_NLP*             nlp,                /**< current NLP data */
   SCIP_VAR*             var                 /**< variable to get solution value for */
   )
{
   int varpos;

   assert(nlp != NULL);
   assert(var != NULL);
   assert(nlp->solstat > SCIP_NLPSOLSTAT_FEASIBLE || nlp->primalsolution != NULL);

   if( nlp->solstat > SCIP_NLPSOLSTAT_FEASIBLE || nlp->primalsolution == NULL )
      return SCIP_INVALID;

   assert(SCIPhashmapExists(nlp->varhash, var));
   varpos = (int) (size_t) SCIPhashmapGetImage(nlp->varhash, var);

   return nlp->primalsolution[varpos];
}

/** gets integer parameter of NLP */
SCIP_RETCODE SCIPnlpGetIntPar(
   SCIP_NLP*             nlp,                /**< pointer to NLP datastructure */
   SCIP_NLPPARAM         type,               /**< parameter number */
   int*                  ival                /**< pointer to store the parameter value */
)
{
   assert(nlp  != NULL);
   assert(nlp->solver  != NULL);
   assert(nlp->problem != NULL);
   assert(ival != NULL);

   SCIP_CALL( SCIPnlpiGetIntPar(nlp->solver, nlp->problem, type, ival) );

   return SCIP_OKAY;
}

/** sets integer parameter of NLP */
SCIP_RETCODE SCIPnlpSetIntPar(
   SCIP_NLP*             nlp,                /**< pointer to NLP datastructure */
   SCIP_NLPPARAM         type,               /**< parameter number */
   int                   ival                /**< parameter value */
)
{
   assert(nlp  != NULL);
   assert(nlp->solver  != NULL);
   assert(nlp->problem != NULL);

   SCIP_CALL( SCIPnlpiSetIntPar(nlp->solver, nlp->problem, type, ival) );

   return SCIP_OKAY;
}

/** gets floating point parameter of NLP */
SCIP_RETCODE SCIPnlpGetRealPar(
   SCIP_NLP*             nlp,                /**< pointer to NLP datastructure */
   SCIP_NLPPARAM         type,               /**< parameter number */
   SCIP_Real*            dval                /**< pointer to store the parameter value */
)
{
   assert(nlp  != NULL);
   assert(nlp->solver  != NULL);
   assert(nlp->problem != NULL);
   assert(dval != NULL);

   SCIP_CALL( SCIPnlpiGetRealPar(nlp->solver, nlp->problem, type, dval) );

   return SCIP_OKAY;
}

/** sets floating point parameter of NLP */
SCIP_RETCODE SCIPnlpSetRealPar(
   SCIP_NLP*             nlp,                /**< pointer to NLP datastructure */
   SCIP_NLPPARAM         type,               /**< parameter number */
   SCIP_Real             dval                /**< parameter value */
)
{
   assert(nlp  != NULL);
   assert(nlp->solver  != NULL);
   assert(nlp->problem != NULL);

   SCIP_CALL( SCIPnlpiSetRealPar(nlp->solver, nlp->problem, type, dval) );

   return SCIP_OKAY;
}

/** gets string parameter of NLP */
SCIP_RETCODE SCIPnlpGetStringPar(
   SCIP_NLP*             nlp,                /**< pointer to NLP datastructure */
   SCIP_NLPPARAM         type,               /**< parameter number */
   const char**          sval                /**< pointer to store the parameter value */
)
{
   assert(nlp  != NULL);
   assert(nlp->solver  != NULL);
   assert(nlp->problem != NULL);
   assert(sval != NULL);

   SCIP_CALL( SCIPnlpiGetStringPar(nlp->solver, nlp->problem, type, sval) );

   return SCIP_OKAY;
}

/** sets string parameter of NLP */
SCIP_RETCODE SCIPnlpSetStringPar(
   SCIP_NLP*             nlp,                /**< pointer to NLP datastructure */
   SCIP_NLPPARAM         type,               /**< parameter number */
   const char*           sval                /**< parameter value */
)
{
   assert(nlp  != NULL);
   assert(nlp->solver  != NULL);
   assert(nlp->problem != NULL);

   SCIP_CALL( SCIPnlpiSetStringPar(nlp->solver, nlp->problem, type, sval) );

   return SCIP_OKAY;
}
