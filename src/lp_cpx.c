/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2002 Tobias Achterberg                              */
/*                            Thorsten Koch                                  */
/*                            Alexander Martin                               */
/*                  2002-2002 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the SCIP Academic Licence.        */
/*                                                                           */
/*  You should have received a copy of the SCIP Academic License             */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   lp_cpx.c
 * @brief  LP interface for CPLEX
 * @author Tobias Achterberg
 */

/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/


#if 0  /* ??? TODO */


#include <stdlib.h>
#include "cplex.h"
#include "lp.h"
#include "def.h"
#include "misc.h"

#define CPX_NOTCALLED  -1

struct cpxlp
{
   CPXLPptr cpxlpptr;    /* pointer to CPLEX LP                    */
   int      cpxretval;   /* return value of last optimization call */
};
typedef struct cpxlp CPXLP;

static int CPX_return_val = CPX_NOTCALLED;

static void
CPX_invalidateSolution( CPXLP *cpxlp )
{
   assert( cpxlp != NULL );
   cpxlp->cpxretval = CPX_NOTCALLED;
}

static int
CPX_isValidSolution( CPXLP *cpxlp )
{
   assert( cpxlp != NULL );
   return( cpxlp->cpxretval == 0
      || cpxlp->cpxretval == CPXERR_PRESLV_INForUNBD
      || cpxlp->cpxretval == CPXERR_PRESLV_INF
      || cpxlp->cpxretval == CPXERR_PRESLV_UNBD);
}

static int
CPX_isInfeasibleSolution( CPXLP *cpxlp )
{
   assert( cpxlp != NULL );
   return( cpxlp->cpxretval == CPXERR_PRESLV_INForUNBD
      || cpxlp->cpxretval == CPXERR_PRESLV_INF );
}

static int
CPX_isUnboundedSolution( CPXLP *cpxlp )
{
   assert( cpxlp != NULL );
   return( cpxlp->cpxretval == CPXERR_PRESLV_INForUNBD
      || cpxlp->cpxretval == CPXERR_PRESLV_UNBD);
}

static int
CPX_isOptimalSolution( CPXLP *cpxlp )
{
   assert( cpxlp != NULL );
   return( cpxlp->cpxretval == 0 );
}

int
SIPstrongbranch( FILE * ferr, SIPInfaLP infaLP, SIPLP lptr, double *psol,
   int *cand, int ncand, double *down, double *up, int itlim )
{
   CPXLP *cpxlp = (CPXLP*)lptr;
   int   restat;

   assert( cpxlp != NULL );
   restat =
      CPXstrongbranch( (CPXENVptr) (void *) infaLP, cpxlp->cpxlpptr,
         cand, ncand, down, up, itlim);
   if( restat != 0 )
   {
      fprintf( ferr, "Error in CPXstrongbranch (), %d returned\n", restat );
      return SIP_LPERROR;
   }

   return SIP_OKAY;
}				/* END STRONGBRANCH */

int
SIPoptLP( SIPInfaLP infaLP, SIPLP lptr )
{
   CPXLP *cpxlp = (CPXLP*)lptr;

   assert( cpxlp != NULL );
   cpxlp->cpxretval =
      CPXdualopt( (CPXENVptr) (void *) infaLP, cpxlp->cpxlpptr );
   switch( cpxlp->cpxretval  )
   {
   case 0:
   case CPXERR_PRESLV_INForUNBD:
   case CPXERR_PRESLV_INF:
   case CPXERR_PRESLV_UNBD:
      return SIP_OKAY;
   case CPXERR_NO_MEMORY:
      return SIP_NOMEMORY;
   default:
      return SIP_LPERROR;
   }
}				/* END OPTLP */

int
SIPoptLPPrimal( SIPInfaLP infaLP, SIPLP lptr )
{
   CPXLP *cpxlp = (CPXLP*)lptr;

   assert( cpxlp != NULL );
   cpxlp->cpxretval =
      CPXprimopt( (CPXENVptr) (void *) infaLP, cpxlp->cpxlpptr );
   switch( cpxlp->cpxretval )
   {
   case 0:
   case CPXERR_PRESLV_INForUNBD:
   case CPXERR_PRESLV_INF:
   case CPXERR_PRESLV_UNBD:
      return SIP_OKAY;
   case CPXERR_NO_MEMORY:
      return SIP_NOMEMORY;
   default:
      return SIP_LPERROR;
   }
}				/* END OPTLPPRIMAL */

int
SIPgetLPStatus( SIPInfaLP infaLP, SIPLP lptr, int *solstat )
{
   CPXLP *cpxlp = (CPXLP*)lptr;

   assert( cpxlp != NULL );
   if( CPX_isValidSolution( cpxlp ) )
   {
      /* first check for primal infeasible (= dual unbounded) solution,
         because a CPXERR_PRESLV_INForUNBD should be considered primal
         infeasible and dual unbounded */
      if( CPX_isUnboundedSolution( cpxlp ) )
	 *solstat = CPX_UNBOUNDED;
      else if( CPX_isInfeasibleSolution( cpxlp ) )
	 *solstat = CPX_INFEASIBLE;
      else if( CPX_isOptimalSolution( cpxlp ) )
      {
	 *solstat =
	    CPXgetstat( (CPXENVptr) (void *) infaLP, cpxlp->cpxlpptr );
	 if( *solstat == 0 )
	    return SIP_LPERROR;
      }
      else
	 return SIP_LPERROR;

      return SIP_OKAY;
   }

   return SIP_LPERROR;
}

/* returns TRUE iff error occured during LP solve */
int
SIPerrorLP( SIPLP lptr, int solstat )
{
   CPXLP *cpxlp = (CPXLP*)lptr;
   return !CPX_isValidSolution( cpxlp );
}				/* END ERRORLP */

/* returns TRUE iff LP solution is stable, ie
   no numerical troubles occured */
int
SIPisStable( SIPLP lptr, int solstat )
{
   CPXLP *cpxlp = (CPXLP*)lptr;
   return( CPX_isValidSolution( cpxlp )
      && solstat != CPX_NUM_BEST_FEAS
      && solstat != CPX_NUM_BEST_INFEAS
      && solstat != CPX_OPTIMAL_INFEAS );
}				/* END ISSTABLE */

/* returns TRUE iff LP is primal infeasible */
int
SIPisPrimalInfeasible( SIPLP lptr, int solstat )
{
   /* primal infeasible means dual unbounded */
   CPXLP *cpxlp = (CPXLP*)lptr;
   return( CPX_isUnboundedSolution( cpxlp )
      || solstat == CPX_UNBOUNDED );
}				/* END ISPRESOLVED */

/* returns TRUE iff LP is primal unbounded */
int
SIPisPrimalUnbounded( SIPLP lptr, int solstat )
{
   /* primal unbounded means dual infeasible */
   CPXLP *cpxlp = (CPXLP*)lptr;
   return( CPX_isInfeasibleSolution( cpxlp )
      || solstat == CPX_INFEASIBLE );
}				/* ISINFEAS */

/* returns TRUE iff LP is solved to optimality */
int
SIPisOptimal( SIPLP lptr, int solstat )
{
   CPXLP *cpxlp = (CPXLP*)lptr;
   return( CPX_isOptimalSolution( cpxlp )
      && solstat == CPX_OPTIMAL );
}				/* ISOPTIMAL */

/* returns TRUE iff a dual feasible solution has been found */
int
SIPisDualValid( SIPLP lptr, int solstat )
{
   CPXLP *cpxlp = (CPXLP*)lptr;
   return( CPX_isValidSolution( cpxlp )
      && ( solstat == CPX_OPTIMAL
	 || solstat == CPX_OBJ_LIM
         || solstat == CPX_IT_LIM_FEAS
	 || solstat == CPX_TIME_LIM_FEAS
         || solstat == CPX_NUM_BEST_FEAS
	 || solstat == CPX_ABORT_FEAS ) );
}				/* ISFEASIBLE */

/* returns TRUE iff objective limit is exceeded */
int
SIPexObjlim( SIPLP lptr, int solstat )
{
   CPXLP *cpxlp = (CPXLP*)lptr;
   return( CPX_isUnboundedSolution( cpxlp )
      || solstat == CPX_UNBOUNDED
      || solstat == CPX_OBJ_LIM );
}				/* END EXOBJLIM */

/* TRUE iff iteration limit has been exceeded */
int
SIPiterlim( SIPLP lptr, int solstat )
{
   return( solstat == CPX_IT_LIM_FEAS
      || solstat == CPX_IT_LIM_INFEAS );
}				/* END ITERLIM */

/* TRUE iff time limit has been exceeded */
int
SIPtimelim( SIPLP lptr, int solstat )
{
   return( solstat == CPX_TIME_LIM_FEAS
      || solstat == CPX_TIME_LIM_INFEAS );
}				/* END TIMELIM */

int
SIPsetbase( FILE * ferr, SIPInfaLP infaLP, SIPLP lptr, double *dnorm,
   int *cstat, int *rstat, int pricing )
{
   CPXENVptr CPXenv = (CPXENVptr) (void *) infaLP;
   CPXLP     *cpxlp = (CPXLP*)lptr;
   int       restat;

   assert( cpxlp != NULL );
   CPX_invalidateSolution( cpxlp );

   if( pricing == CPX_DPRIIND_STEEP && dnorm != NULL )
   {
      restat = CPXcopybasednorms( CPXenv, cpxlp->cpxlpptr, cstat, rstat, dnorm );
      if( restat != 0 )
      {
	 fprintf(ferr, "Error in CPXloaddnorms (), %d returned\n", restat);
	 return SIP_LPERROR;
      }
   }
   else
   {
      restat = CPXcopybase( CPXenv, cpxlp->cpxlpptr, cstat, rstat );
      if( restat != 0 )
      {
	 fprintf(ferr, "Error in CPXloadbase (), %d returned\n", restat);
	 return SIP_LPERROR;
      }
   }

   return SIP_OKAY;
}				/* END SETBASE */

int
SIPgetsol(FILE * ferr, SIPInfaLP infaLP, SIPLP lptr, double *objval,
   double *psol, double *pi, double *slck, double *redcost)
{
   CPXLP   *cpxlp = (CPXLP*)lptr;
   int     restat;
   int     dummy;

   assert( cpxlp != NULL );

   restat =
      CPXsolution( (CPXENVptr) (void *) infaLP, cpxlp->cpxlpptr,
         &dummy, objval, psol, pi, slck, redcost );
   if( restat != 0 )
   {
      fprintf(ferr, "Error in CPXsolution (), %d returned.\n", restat);
      return SIP_LPERROR;
   }
   
   return SIP_OKAY;
}				/* END GETSOL */

int
SIPgetbase(FILE * ferr, SIPInfaLP infaLP, SIPLP lptr, double *dnorm,
   int *cstat, int *rstat)
{
   CPXLP   *cpxlp = (CPXLP*)lptr;
   int     restat;

   assert( cpxlp != NULL );
   if( dnorm != NULL )
   {
      restat =
	 CPXgetbasednorms( (CPXENVptr) (void *) infaLP,
            cpxlp->cpxlpptr, cstat, rstat, dnorm );
      if( restat != 0 )
      {
	 fprintf(ferr, "Error in CPXgetbasednorms (), returned %d\n", restat);
	 return SIP_LPERROR;
      }
   }
   else
   {
      restat =
	 CPXgetbase( (CPXENVptr) (void *) infaLP, cpxlp->cpxlpptr,
            cstat, rstat );
      if( restat != 0 )
      {
	 fprintf(ferr, "Error in CPXgetbase (), returned %d\n", restat);
	 return SIP_LPERROR;
      }
   }

   return SIP_OKAY;
}				/* END GETBASE */

void
SIPsetintparLP( SIPInfaLP infaLP, SIPLP lptr, int type, int ival )
{
   /* WARNING! CPLEX sets parameters in its environment, which
    * affects all LPs
    */

   CPXENVptr CPXenv = (CPXENVptr) (void *) infaLP;
   CPXLP     *cpxlp = (CPXLP*)lptr;

   assert( infaLP != NULL );
   assert( cpxlp != NULL );

   switch( type )
   {
   case SIP_FROMSCRATCH:
      CPXsetintparam( CPXenv, CPX_PARAM_ADVIND,
	 (ival == SIP_OFF) ? CPX_ON : CPX_OFF );
      break;
   case SIP_LPITLIM:
      CPXsetintparam( CPXenv, CPX_PARAM_ITLIM, ival );
      break;
   case SIP_FASTMIP:
      assert( (ival == SIP_ON) || (ival == SIP_OFF) );
      CPXsetintparam( CPXenv, CPX_PARAM_FASTMIP,
	 (ival == SIP_ON) ? CPX_ON : CPX_OFF );
      break;
   case SIP_PRICING:
      switch( ival )
      {
      case SIP_FULLPRICE:
	 CPXsetintparam( CPXenv, CPX_PARAM_DPRIIND, CPX_DPRIIND_FULL );
	 break;
      case SIP_STEEPEDGE:
	 CPXsetintparam( CPXenv, CPX_PARAM_DPRIIND, CPX_DPRIIND_STEEP );
	 break;
      case SIP_APPROXSTE:
	 CPXsetintparam( CPXenv, CPX_PARAM_DPRIIND, CPX_DPRIIND_STEEPQSTART );
	 break;
      default:
         abort;
      }
      break;
   case SIP_LPINFO:
      if( ival == SIP_ON )
      {
	 CPXsetintparam( CPXenv, CPX_PARAM_SIMDISPLAY, CPX_ON );
	 CPXsetintparam( CPXenv, CPX_PARAM_SCRIND, CPX_ON );
      }
      else 
      {
         assert(ival == SIP_OFF);
         
	 CPXsetintparam( CPXenv, CPX_PARAM_SIMDISPLAY, CPX_OFF );
	 CPXsetintparam( CPXenv, CPX_PARAM_SCRIND, CPX_OFF );
      }
      break;
   default:
      abort();
   }
}

void
SIPsetdblparLP( SIPInfaLP infaLP, SIPLP lptr, int type, double dval )
{
   /* WARNING! CPLEX sets parameters in its environment, which
    * affects all LPs
    */

   CPXENVptr CPXenv = (CPXENVptr) (void *) infaLP;
   CPXLP     *cpxlp = (CPXLP*)lptr;

   assert( infaLP != NULL );
   assert( cpxlp != NULL );

   switch( type )
   {
   case SIP_FEASTOL:
      CPXsetdblparam( CPXenv, CPX_PARAM_EPRHS, dval );
      break;
   case SIP_LOBJLIM:
      CPXsetdblparam( CPXenv, CPX_PARAM_OBJLLIM, dval );
      break;
   case SIP_UOBJLIM:
      CPXsetdblparam( CPXenv, CPX_PARAM_OBJULIM, dval );
      break;
   case SIP_LPTILIM:
      CPXsetdblparam( CPXenv, CPX_PARAM_TILIM, dval );
      break;
   default:
      abort();
   }
}

int
SIPgetintparLP( SIPInfaLP infaLP, SIPLP lptr, int type, int *ival )
{
   CPXENVptr CPXenv = (CPXENVptr) (void *) infaLP;
   CPXLP     *cpxlp = (CPXLP*)lptr;
   int       restat = SIP_OKAY;
   int       advind;

   assert( infaLP != NULL );
   assert( cpxlp != NULL );
   assert( ival != NULL );

   switch( type )
   {
   case SIP_FROMSCRATCH:
      CPXgetintparam( CPXenv, CPX_PARAM_ADVIND, &advind );
      if( advind == CPX_ON )
	 *ival = SIP_OFF;
      else
	 *ival = SIP_ON;
      break;
   case SIP_LPNROW:
      assert( cpxlp->cpxlpptr != NULL );
      *ival = CPXgetnumrows( CPXenv, cpxlp->cpxlpptr );
      break;
   case SIP_LPIT:
      assert( cpxlp->cpxlpptr != NULL );
      *ival = CPXgetitcnt( CPXenv, cpxlp->cpxlpptr );
      break;
   case SIP_LPIT1:
      assert( cpxlp->cpxlpptr != NULL );
      *ival = CPXgetphase1cnt( CPXenv, cpxlp->cpxlpptr );
      break;
   default:
      restat = SIP_LPERROR;
      break;
   }
   return restat;
}				/* END GETPARLP */

int
SIPgetdblparLP( SIPInfaLP infaLP, SIPLP lptr, int type, double *dval )
{
   CPXENVptr CPXenv = (CPXENVptr) (void *) infaLP;
   CPXLP     *cpxlp = (CPXLP*)lptr;
   int       restat = SIP_OKAY;
   int       status;

   assert( infaLP != NULL );
   assert( cpxlp != NULL );
   assert( dval != NULL );

   switch( type )
   {
   case SIP_FEASTOL:
      status = CPXgetdblparam( CPXenv, CPX_PARAM_EPRHS, dval );
      if( status != 0 )
	 restat = SIP_LPERROR;
      break;
   default:
      restat = SIP_LPERROR;
      break;
   }
   return restat;
}				/* END GETPARLP */

int
SIPgetDefaultLPFeastol( SIPInfaLP infaLP, double* LPFeastol )
{
   CPXENVptr CPXenv = (CPXENVptr) (void *) infaLP;
   int       restat = SIP_OKAY;
   int       status;

   status = CPXgetdblparam( CPXenv, CPX_PARAM_EPRHS, LPFeastol );
   if( status != 0 )
      restat = SIP_LPERROR;
   return restat;
}

int
SIPwriteLP( FILE * ferr, SIPInfaLP infaLP, SIPLP lptr, char *fname )
{
   CPXLP   *cpxlp = (CPXLP*)lptr;
   int     restat;

   assert( cpxlp != NULL );
   restat =
      CPXlpwrite( (CPXENVptr) (void *) infaLP, cpxlp->cpxlpptr, fname );
   if( restat != 0 )
   {
      fprintf( ferr, "Error in CPXlpwrite, %d returned\n", restat );
      return SIP_LPERROR;
   }

   return SIP_OKAY;
}				/* END WRITELP */

int
SIPwriteB( FILE * ferr, SIPInfaLP infaLP, SIPLP lptr, char *fname )
{
   CPXLP   *cpxlp = (CPXLP*)lptr;
   int     restat;

   assert( cpxlp != NULL );
   restat =
      CPXmbasewrite( (CPXENVptr) (void *) infaLP, cpxlp->cpxlpptr, fname );
   if( restat != 0 )
   {
      fprintf( ferr, "Error in CPXmbasewrite (), %d returned\n", restat );
      return SIP_LPERROR;
   }

   return SIP_OKAY;
}				/* END WRITEB */


int
SIPgetlb( FILE * ferr, SIPInfaLP infaLP, SIPLP lptr, double *lb, int beg,
   int end )
{
   CPXLP   *cpxlp = (CPXLP*)lptr;
   int     restat;

   assert( cpxlp != NULL );
   restat =
      CPXgetlb( (CPXENVptr) (void *) infaLP, cpxlp->cpxlpptr, lb, beg, end );
   if( restat != 0 )
   {
      fprintf( ferr, "Error in CPXgetlb (), %d returned\n", restat );
      return SIP_LPERROR;
   }

   return SIP_OKAY;
}				/* END GETLB */

int
SIPgetub( FILE * ferr, SIPInfaLP infaLP, SIPLP lptr, double *ub, int beg,
   int end )
{
   CPXLP   *cpxlp = (CPXLP*)lptr;
   int     restat;

   assert( cpxlp != NULL );
   restat =
      CPXgetub( (CPXENVptr) (void *) infaLP, cpxlp->cpxlpptr, ub, beg, end );
   if( restat != 0 )
   {
      fprintf( ferr, "Error in CPXgetlb (), %d returned\n", restat );
      return SIP_LPERROR;
   }

   return SIP_OKAY;
}				/* END GETUB */

int
SIPchgbds( FILE * ferr, SIPInfaLP infaLP, SIPLP lptr, int n, int *ind, char *lu,
   double *bd )
{
   CPXLP   *cpxlp = (CPXLP*)lptr;
   int     restat;

   assert( cpxlp != NULL );

   CPX_invalidateSolution( cpxlp );

   restat =
      CPXchgbds( (CPXENVptr) (void *) infaLP, cpxlp->cpxlpptr, n, ind, lu, bd );
   if( restat != 0 )
   {
      fprintf( ferr, "Error in CPXchgbds (), %d returned\n", restat );
      return SIP_LPERROR;
   }

   return SIP_OKAY;
}				/* END CHGBDS */

int
SIPchgrhs( FILE * ferr, SIPInfaLP infaLP, SIPLP lptr, int n, int *ind, 
   double *rhs )
{
   CPXLP   *cpxlp = (CPXLP*)lptr;
   int     restat;

   assert( cpxlp != NULL );

   CPX_invalidateSolution( cpxlp );

   restat =
      CPXchgrhs( (CPXENVptr) (void *) infaLP, cpxlp->cpxlpptr, n, ind, rhs );
   if( restat != 0 )
   {
      fprintf( ferr, "Error in CPXchgrhs (), %d returned\n", restat );
      return SIP_LPERROR;
   }

   return SIP_OKAY;
}				/* END CHGRHS */

void
SIPchgobjsen( FILE * ferr, SIPInfaLP infaLP, SIPLP lptr, int objsen )
{
   CPXLP   *cpxlp = (CPXLP*)lptr;

   assert( cpxlp != NULL );

   CPX_invalidateSolution( cpxlp );

   CPXchgobjsen( (CPXENVptr) (void *) infaLP, cpxlp->cpxlpptr, objsen );
}				/* END CHGOBJSEN */

int
SIPdelrows( FILE * ferr, SIPInfaLP infaLP, SIPLP lptr, int *dstat )
{
   CPXLP   *cpxlp = (CPXLP*)lptr;
   int     restat;

   assert( cpxlp != NULL );

   CPX_invalidateSolution( cpxlp );

   restat =
      CPXdelsetrows( (CPXENVptr) (void *) infaLP, cpxlp->cpxlpptr,
         dstat );
   if( restat != 0 )
   {
      fprintf(ferr, "Error in CPXdelsetrows (), %d returned\n", restat);
      return SIP_LPERROR;
   }

   return SIP_OKAY;
}				/* END DELROWS */

int
SIPgetBind( FILE * ferr, SIPInfaLP infaLP, SIPLP lptr, int *head )
{
   CPXLP   *cpxlp = (CPXLP*)lptr;
   int     restat;

   assert( cpxlp != NULL );

   restat =
      CPXgetbhead( (CPXENVptr) (void *) infaLP, cpxlp->cpxlpptr, head,
         NULL );
   if( restat != 0 )
   {
      fprintf( ferr, "Error in CPXgetbhead(), %d returned.\n", restat );
      return SIP_LPERROR;
   }

   return SIP_OKAY;
}				/* END GETBIND */

int
SIPgetrow( FILE * ferr, SIPInfaLP infaLP, SIPLP lptr, int i, double *val,
   int *ind, int *nnonz )
{
   CPXLP   *cpxlp = (CPXLP*)lptr;
   int     restat;
   int     beg[1];
   int     surplus;

   assert( cpxlp != NULL );
   restat =
      CPXgetrows( (CPXENVptr) (void *) infaLP, cpxlp->cpxlpptr, nnonz,
         beg, ind, val, *nnonz, &surplus, i, i );
   if( restat != 0 )
   {
      fprintf( ferr, "Error in CPXgetrows(), %d returned.\n", restat );
      return SIP_LPERROR;
   }
   *nnonz -= surplus;

   return SIP_OKAY;
}				/* END GETROW */

int
SIPgetrowBinv( FILE * ferr, SIPInfaLP infaLP, SIPLP lptr, int i, double *val )
{
   CPXLP   *cpxlp = (CPXLP*)lptr;
   int     restat;

   assert( cpxlp != NULL );
   restat =
      CPXbinvrow( (CPXENVptr) (void *) infaLP, cpxlp->cpxlpptr, i,
         val );
   if( restat != 0 )
   {
      fprintf( ferr, "Error in CPXbinvrow(), %d returned.\n", restat );
      return SIP_LPERROR;
   }

   return SIP_OKAY;
}				/* END GETROWBINV */

int
SIPgetrowBinvA( FILE * ferr, SIPInfaLP infaLP, SIPLP lptr, int i, double *binv,
   double *val )
{
   CPXLP   *cpxlp = (CPXLP*)lptr;
   int     restat;

   assert( cpxlp != NULL );
   restat =
      CPXbinvarow( (CPXENVptr) (void *) infaLP, cpxlp->cpxlpptr, i,
         val );
   if( restat != 0 )
   {
      fprintf( ferr, "Error in CPXbinarow(), %d returned.\n", restat );
      return SIP_LPERROR;
   }

   return SIP_OKAY;
}				/* END GETROWBINVA */

int
SIPcopyLP( FILE * ferr, SIPInfaLP infaLP, SIPLP lptr, int ncol, int nrow,
   int objsen, double *obj, double *rhs, char *sen, int *beg, int *cnt,
   int *ind, double *val, double *lb, double *ub, char **cname, char **rname )
{
   CPXLP   *cpxlp = (CPXLP*)lptr;
   int     restat;

   assert( cpxlp != NULL );
   restat =
      CPXcopylpwnames( (CPXENVptr) (void *) infaLP, cpxlp->cpxlpptr,
         ncol, nrow, objsen, obj, rhs, sen, beg, cnt, ind, val, lb, ub, NULL,
         cname, rname );
   if( restat != 0 )
   {
      fprintf( ferr, "Copying LP failed, %d returned\n", restat );
      return SIP_LPERROR;
   }

   return SIP_OKAY;
}				/* END COPYLP */

int
SIPopenLP( FILE * ferr, SIPInfaLP infaLP, SIPLP *lptr, char *name )
{
   CPXLP   *cpxlp;
   int     restat;

   *lptr = NULL;

   ALLOC_OKAY( allocMemory(cpxlp) );

   CPX_invalidateSolution( cpxlp );

   cpxlp->cpxlpptr = CPXcreateprob( (CPXENVptr) (void *) infaLP, &restat, name );
   if( restat != 0 )
   {
      fprintf( ferr, "Creating problem failed, %d returned\n", restat );
      freeMemory(cpxlp);
      return SIP_LPERROR;
   }

   *lptr = (SIPLP)cpxlp;
   return SIP_OKAY;
}				/* END OPENLP */

int
SIPfreeLP( FILE * ferr, SIPInfaLP infaLP, SIPLP *lptr )
{
   CPXLP **cpxlp = (CPXLP**)lptr;

   if( *cpxlp != NULL )
   {
      int     restat;

      restat = CPXfreeprob( (CPXENVptr) (void *) infaLP, &((*cpxlp)->cpxlpptr) );
      if( restat != 0 )
      {
	 fprintf( ferr, "CPXfreeprob () failed, %d returned\n", restat );
	 return SIP_LPERROR;
      }
      freeMemory(*cpxlp);
      *cpxlp = NULL;
   }
   return SIP_OKAY;
}				/* END FREELP */

int
SIPfreeInfa( FILE * ferr, SIPInfaLP *infaLP )
{
   int     restat;
   
   /* close LP interface */
   if( *infaLP != NULL )
   {
      restat = CPXcloseCPLEX( (CPXENVptr *) (void **) infaLP );
      if( restat != 0 )
      {
         fprintf( ferr, "Error in CPXcloseCPLEX (), %d returned\n", restat );
         return SIP_LPERROR;
      }
   }
   return SIP_OKAY;
}				/* END FREEINFA */

int
SIPopenInfa( FILE * ferr, SIPInfaLP *infaLP )
{
   int restat;

   *infaLP = (SIPInfaLP) (void *) CPXopenCPLEX( &restat );
   if( restat != 0 )
   {
      fprintf( ferr, "Error in CPXopenCPLEX (), %d returned.\n", restat );
      return SIP_LPERROR;
   }
   return SIP_OKAY;
}				/* OPENINFA */

int
SIPaddrow( FILE * ferr, SIPInfaLP infaLP, SIPLP lptr, int nrow, int nnonz,
   double *rhs, char *sen, int *beg, int *ind, double *val, char **name )
{
   CPXLP   *cpxlp = (CPXLP*)lptr;
   int     restat;

   assert( cpxlp != NULL );

   CPX_invalidateSolution( cpxlp);

   restat =
      CPXaddrows( (CPXENVptr) (void *) infaLP, cpxlp->cpxlpptr, 0,
         nrow, nnonz, rhs, sen, beg, ind, val, NULL, name );
   if( restat != 0 )
   {
      fprintf( ferr, "Error in CPXaddrows (), %d returned\n", restat );
      return SIP_LPERROR;
   }

   return SIP_OKAY;
}				/* END ADDROW */

int
SIPreadFile( SET *set, MIP *mip, const char *filename )
{
   CPXENVptr cpxenv = NULL;
   CPXLPptr cpxlpptr = NULL;
   SIPLP   siplp = NULL;
   int     restat;
   int     sipstat;
   int     i;
   int     numnz;
   int     retnumnz;
   int     surplus;

   assert(set != NULL && mip != NULL && filename != NULL);

   restat = 0;
   sipstat = SIP_LPERROR;

   /* check filename */
   if( strlen(filename) == 0 )
   {
      fprintf(stderr, "No file specified\n");
      sipstat = SIP_NOFILE;
      goto ERROR;
   }
   /* open CPLEX LP (only for reading data) */
   restat = SIPopenLP(stderr, set->inface->lp, &siplp, (char*)filename);
   if( restat != 0 )
   {
      sipstat = SIP_LPERROR;
      goto ERROR;
   }
   cpxenv = (CPXENVptr) (void *) (set->inface->lp);
   cpxlpptr = ((CPXLP*)siplp)->cpxlpptr;

   /* store filename in MIP structure */
   allocMemoryArray(mip->name, strlen(filename)+1);
   if( mip->name == NULL )
   {
      sipstat = SIP_NOMEMORY;
      goto ERROR;
   }
   strcpy(mip->name, filename);

   /* read file into CPLEX */
   restat = CPXreadcopyprob(cpxenv, cpxlpptr, (char*)filename, NULL);
   if( restat != 0 )
   {
      sipstat = SIP_NOFILE;
      goto ERROR;
   }

   /* get matrix size */
   numnz = CPXgetnumnz(cpxenv, cpxlpptr);

   /* get number of columns */
   mip->ncol = CPXgetnumcols(cpxenv, cpxlpptr);

   /* get matrix in column form */
   allocMemory(mip->col);
   if( mip->col == NULL )
   {
      sipstat = SIP_NOMEMORY;
      goto ERROR;
   }
   allocMemoryArray(mip->col->pnt, mip->ncol + 1);
   allocMemoryArray(mip->col->ind, numnz);
   allocMemoryArray(mip->col->val, numnz);
   if( mip->col->pnt == NULL || mip->col->ind == NULL || mip->col->val == NULL )
   {
      sipstat = SIP_NOMEMORY;
      goto ERROR;
   }

   restat = CPXgetcols(cpxenv, cpxlpptr, &retnumnz, mip->col->pnt, mip->col->ind,
      mip->col->val, numnz, &surplus, 0, mip->ncol - 1);
   if( restat != 0 )
      goto ERROR;
   assert(retnumnz == numnz && surplus == 0);
   mip->col->pnt[mip->ncol] = numnz;

   /* get objective function */
   allocMemoryArray(mip->obj, mip->ncol);
   if( mip->obj == NULL )
   {
      sipstat = SIP_NOMEMORY;
      goto ERROR;
   }

   restat = CPXgetobj(cpxenv, cpxlpptr, mip->obj, 0, mip->ncol - 1);
   if( restat != 0 )
      goto ERROR;

   /* get lower bounds */
   allocMemoryArray(mip->lb, mip->ncol);
   if( mip->lb == NULL )
   {
      sipstat = SIP_NOMEMORY;
      goto ERROR;
   }

   restat = CPXgetlb(cpxenv, cpxlpptr, mip->lb, 0, mip->ncol - 1);
   if( restat != 0 )
      goto ERROR;

   /* get upper bounds */
   allocMemoryArray(mip->ub, mip->ncol);
   if( mip->ub == NULL )
   {
      sipstat = SIP_NOMEMORY;
      goto ERROR;
   }

   restat = CPXgetub(cpxenv, cpxlpptr, mip->ub, 0, mip->ncol - 1);
   if( restat != 0 )
      goto ERROR;

   /* get type of variables */
   allocMemoryArray(mip->ctype, mip->ncol);
   if( mip->ctype == NULL )
   {
      sipstat = SIP_NOMEMORY;
      goto ERROR;
   }

   restat = CPXgetctype(cpxenv, cpxlpptr, mip->ctype, 0, mip->ncol - 1);
   if( restat != 0 )
      goto ERROR;

   /* get variable names */
   allocMemoryArray(mip->cname, mip->ncol);
   if( mip->cname == NULL )
   {
      sipstat = SIP_NOMEMORY;
      goto ERROR;
   }

   restat = CPXgetcolname(cpxenv, cpxlpptr, mip->cname, NULL, 0, &surplus, 0, mip->ncol - 1);
   if( restat != 0 && restat != CPXERR_NEGATIVE_SURPLUS )
      goto ERROR;

   assert(surplus <= 0);
   mip->cstoresz = -surplus;
   allocMemoryArray(mip->cstore, mip->cstoresz);
   if( mip->cstore == NULL )
   {
      sipstat = SIP_NOMEMORY;
      goto ERROR;
   }

   restat = CPXgetcolname(cpxenv, cpxlpptr, mip->cname, mip->cstore, mip->cstoresz, &surplus, 0, mip->ncol - 1);
   if( restat != 0 )
      goto ERROR;
   assert(surplus == 0);

   /* get number of rows */
   mip->nrow = CPXgetnumrows(cpxenv, cpxlpptr);

#if 0	    /* The unpresolved matrix is not needed in row form */
   /* get matrix in row form */
   allocMemory(mip->row);
   if( mip->row == NULL )
   {
      sipstat = SIP_NOMEMORY;
      goto ERROR;
   }
   allocMemoryArray(mip->row->pnt, mip->nrow + 1);
   allocMemoryArray(mip->row->ind, numnz);
   allocMemoryArray(mip->row->val, numnz);
   if( mip->row->pnt == NULL || mip->row->ind == NULL || mip->row->val == NULL )
   {
      sipstat = SIP_NOMEMORY;
      goto ERROR;
   }
   restat = CPXgetrows(cpxenv, cpxlpptr, &retnumnz, mip->row->pnt, mip->row->ind,
      mip->row->val, numnz, &surplus, 0, mip->nrow - 1);
   if( restat != 0 )
      goto ERROR;
   assert(retnumnz == numnz && surplus == 0);
   mip->row->pnt[mip->nrow] = numnz;
#else
   mip->row = NULL;
#endif

   /* get objective sense */
   mip->objsen = CPXgetobjsen(cpxenv, cpxlpptr);
   mip->objoff = 0.0;

   /* get right hand sides */
   allocMemoryArray(mip->rhs, mip->nrow);
   if( mip->rhs == NULL )
   {
      sipstat = SIP_NOMEMORY;
      goto ERROR;
   }

   restat = CPXgetrhs(cpxenv, cpxlpptr, mip->rhs, 0, mip->nrow - 1);
   if( restat != 0 )
      goto ERROR;

   /* get constraint senses */
   allocMemoryArray(mip->sen, mip->nrow);
   if( mip->sen == NULL )
   {
      sipstat = SIP_NOMEMORY;
      goto ERROR;
   }

   restat = CPXgetsense(cpxenv, cpxlpptr, mip->sen, 0, mip->nrow - 1);
   if( restat != 0 )
      goto ERROR;

   /* get row names */
   allocMemoryArray(mip->rname, mip->nrow);
   if( mip->rname == NULL )
   {
      sipstat = SIP_NOMEMORY;
      goto ERROR;
   }

   restat = CPXgetrowname(cpxenv, cpxlpptr, mip->rname, NULL, 0, &surplus, 0, mip->nrow - 1);
   if( restat != 0 && restat != CPXERR_NEGATIVE_SURPLUS )
      goto ERROR;

   assert(surplus <= 0);
   mip->rstoresz = -surplus;
   allocMemoryArray(mip->rstore, mip->rstoresz);
   if( mip->rstore == NULL )
   {
      sipstat = SIP_NOMEMORY;
      goto ERROR;
   }

   restat = CPXgetrowname(cpxenv, cpxlpptr, mip->rname, mip->rstore, mip->rstoresz, &surplus, 0, mip->nrow - 1);
   if( restat != 0 )
      goto ERROR;
   assert(surplus == 0);

   /* get range values for ranged rows */
   allocMemoryArray(mip->rngval, mip->nrow);
   if( mip->rngval == NULL )
   {
      sipstat = SIP_NOMEMORY;
      goto ERROR;
   }

   restat = CPXgetrngval(cpxenv, cpxlpptr, mip->rngval, 0, mip->nrow - 1);
   if( restat == CPXERR_NO_RNGVAL )
      for( i = 0; i < mip->nrow; ++i )
	 mip->rngval[i] = 0.0;
   else if( restat != 0 )
      goto ERROR;

   SIPfreeLP(stderr, set->inface->lp, &siplp);
   return SIP_OKAY;

 ERROR:
   if( restat != 0 )
      printf( "Error (CPLEX:%d) reading file <%s>!\n", restat, filename );
   else
      printf( "Error (SIP:%d) reading file <%s>!\n", sipstat, filename );
   if( siplp != NULL )
      SIPfreeLP(stderr, set->inface->lp, &siplp);

   return sipstat;
}


#endif
