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

#include <assert.h>

#include "cplex.h"
#include "bitencode.h"
#include "lpi.h"


#define CHECK_ZERO(x) { int _restat_; if( (_restat_ = (x)) != 0 ) { errorMessage("LP Error"); \
                                                                    printf("-> CPLEX returned %d.\n", _restat_); \
                                                                    return SCIP_LPERROR; } }
#define NOTCALLED  -1


typedef DUALPACKET COLPACKET;           /* each column needs two bits of information (basic/on_lower/on_upper) */
#define COLS_PER_PACKET DUALPACKETSIZE
typedef SINGLEPACKET ROWPACKET;         /* each row needs one bit of information (basic/nonbasic) */
#define ROWS_PER_PACKET SINGLEPACKETSIZE

/* CPLEX parameter lists which can be changed */
#define NUMINTPARAM  6
static const int    intparam[NUMINTPARAM] = {
   CPX_PARAM_ADVIND,
   CPX_PARAM_ITLIM,
   CPX_PARAM_FASTMIP,
   CPX_PARAM_DPRIIND,
   CPX_PARAM_SIMDISPLAY,
   CPX_PARAM_SCRIND
};
#define NUMDBLPARAM  1
static const double dblparam[NUMDBLPARAM] = {
   CPX_PARAM_EPRHS
};

struct CPXParam                         /**< CPLEX parameter settings */
{
   int              intparval[NUMINTPARAM]; /**< integer parameter values */
   double           dblparval[NUMDBLPARAM]; /**< double parameter values */
};
typedef struct CPXParam CPXPARAM;

enum OptAlgo                            /**< LP algorithm type */
{
   INVALID       = 0,                   /**< problem is not optimized */
   PRIMALSIMPLEX = 1,                   /**< primal simplex algorithm */
   DUALSIMPLEX   = 2                    /**< dual simplex algorithm */
};
typedef enum OptAlgo OPTALGO;

struct LPi                              /**< LP interface */
{
   CPXLPptr         cpxlp;              /**< CPLEX LP pointer */
   OPTALGO          optalgo;            /**< last optimization algorithm */
   int              retval;             /**< return value of last optimization call */
   int              solstat;            /**< solution status of last optimization call */
   CPXPARAM         cpxparam;           /**< actual parameter values for this LP */
};

struct LPState                          /**< LP state stores basis information */
{
   int              numuses:24;         /**< number of times, this LP state is referenced */
   unsigned int     ncol:20;            /**< number of LP columns */
   unsigned int     nrow:20;            /**< number of LP rows */
   COLPACKET*       packcstat;          /**< column basis status in compressed form */
   ROWPACKET*       packrstat;          /**< row basis status in compressed form */
   double*          dnorm;              /**< dual norms of variables */
};


static CPXENVptr    cpxenv = NULL;      /**< CPLEX environment */
static CPXPARAM     defparam;           /**< default CPLEX parameters */
static CPXPARAM     actparam;           /**< actual CPLEX parameters in the environment */
static int          numlp = 0;          /**< number of open LP objects */



static
RETCODE getParameterValues(CPXPARAM* cpxparam)
{
   int i;
   
   assert(cpxenv != NULL);
   assert(cpxparam != NULL);

   for( i = 0; i < NUMINTPARAM; ++i )
   {
      CHECK_ZERO( CPXgetintparam(cpxenv, intparam[i], &(cpxparam->intparval[i])) );
   }
   for( i = 0; i < NUMDBLPARAM; ++i )
   {
      CHECK_ZERO( CPXgetdblparam(cpxenv, dblparam[i], &(cpxparam->dblparval[i])) );
   }

   return SCIP_OKAY;
}
   
static
void checkParameterValues(void)
{
#ifndef NDEBUG
   CPXPARAM par;
   int i;
   
   getParameterValues(&par);
   for( i = 0; i < NUMINTPARAM; ++i )
      assert(actparam.intparval[i] == par.intparval[i]);
   for( i = 0; i < NUMDBLPARAM; ++i )
      assert(actparam.dblparval[i] == par.dblparval[i]);
#endif
}

static
RETCODE setParameterValues(const CPXPARAM* cpxparam)
{
   int i;
   
   assert(cpxenv != NULL);
   assert(cpxparam != NULL);
   
   for( i = 0; i < NUMINTPARAM; ++i )
   {
      if( actparam.intparval[i] != cpxparam->intparval[i] )
      {
         actparam.intparval[i] = cpxparam->intparval[i];
         CHECK_ZERO( CPXsetintparam(cpxenv, intparam[i], actparam.intparval[i]) );
      }
   }
   for( i = 0; i < NUMDBLPARAM; ++i )
   {
      if( actparam.dblparval[i] != cpxparam->dblparval[i] )
      {
         actparam.dblparval[i] = cpxparam->dblparval[i];
         CHECK_ZERO( CPXsetdblparam(cpxenv, dblparam[i], actparam.dblparval[i]) );
      }
   }

   checkParameterValues();

   return SCIP_OKAY;
}

static
void copyParameterValues(CPXPARAM* dest, const CPXPARAM* source)
{
   int i;

   for( i = 0; i < NUMINTPARAM; ++i )
      dest->intparval[i] = source->intparval[i];
   for( i = 0; i < NUMDBLPARAM; ++i )
      dest->dblparval[i] = source->dblparval[i];
}

static
int getIntParam(LPI* lpi, const int param)
{
   int i;
   
   assert(lpi != NULL);

   for( i = 0; i < NUMINTPARAM; ++i )
      if( intparam[i] == param )
         return lpi->cpxparam.intparval[i];

   errorMessage("Unknown CPLEX integer parameter");
   abort();
}

static
double getDblParam(LPI* lpi, const int param)
{
   int i;

   assert(lpi != NULL);

   for( i = 0; i < NUMDBLPARAM; ++i )
      if( dblparam[i] == param )
         return lpi->cpxparam.dblparval[i];

   errorMessage("Unknown CPLEX double parameter");
   abort();
}

static
void setIntParam(LPI* lpi, const int param, int parval)
{
   int i;

   assert(lpi != NULL);

   for( i = 0; i < NUMINTPARAM; ++i )
      if( intparam[i] == param )
      {
         lpi->cpxparam.intparval[i] = parval;
         return;
      }

   errorMessage("Unknown CPLEX integer parameter");
   abort();
}

static
void setDblParam(LPI* lpi, const int param, double parval)
{
   int i;

   assert(lpi != NULL);

   for( i = 0; i < NUMDBLPARAM; ++i )
      if( dblparam[i] == param )
      {
         lpi->cpxparam.dblparval[i] = parval;
         return;
      }

   errorMessage("Unknown CPLEX double parameter");
   abort();
}

static
void invalidateSolution(LPI* lpi)
{
   assert(lpi != NULL);
   lpi->optalgo = INVALID;
   lpi->retval = NOTCALLED;
   lpi->solstat = -1;
}

static int
isValidSolution(LPI* lpi)
{
   assert(lpi != NULL);
   return( lpi->retval == 0
      || lpi->retval == CPXERR_PRESLV_INForUNBD
      || lpi->retval == CPXERR_PRESLV_INF
      || lpi->retval == CPXERR_PRESLV_UNBD);
}

static int
isInfeasibleSolution(LPI* lpi)
{
   assert(lpi != NULL);
   return( lpi->retval == CPXERR_PRESLV_INForUNBD
      || lpi->retval == CPXERR_PRESLV_INF );
}

static int
isUnboundedSolution(LPI* lpi)
{
   assert(lpi != NULL);
   return( lpi->retval == CPXERR_PRESLV_INForUNBD
      || lpi->retval == CPXERR_PRESLV_UNBD);
}

static int
isOptimalSolution(LPI* lpi)
{
   assert(lpi != NULL);
   return( lpi->retval == 0 );
}

static
int cpxObjsen(OBJSEN objsen)
{
   switch( objsen )
   {
   case SCIP_OBJSEN_MAXIMIZE:
      return CPX_MIN;
   case SCIP_OBJSEN_MINIMIZE:
      return CPX_MAX;
   default:
      errorMessage("invalid objective sense");
      abort();
   }
}



RETCODE SCIPopenLP(                     /**< creates an LP problem object */
   LPI**            lpi,                /**< pointer to an LP interface structure */
   const char*      name                /**< name of the LP */
   )
{
   int     restat;

   assert(lpi != NULL);
   assert(numlp >= 0);

   /* create environment */
   if( cpxenv == NULL )
   {
      assert(numlp == 0);
      cpxenv = CPXopenCPLEX(&restat);
      CHECK_ZERO(restat);

      /* get default parameter values */
      getParameterValues(&defparam);
      copyParameterValues(&actparam, &defparam);
   }
   assert(cpxenv != NULL);

   /* create LP */
   allocMemory(*lpi);
   (*lpi)->cpxlp = CPXcreateprob(cpxenv, &restat, (char*)name);
   CHECK_ZERO(restat);
   invalidateSolution(*lpi);
   copyParameterValues(&((*lpi)->cpxparam), &defparam);
   numlp++;

   return SCIP_OKAY;
}

RETCODE SCIPcloseLP(                    /**< deletes an LP problem object */
   LPI**            lpi                 /**< pointer to an LP interface structure */
   )
{
   assert(cpxenv != NULL);
   assert(lpi != NULL);
   assert(*lpi != NULL);

   /* free LP */
   CHECK_ZERO( CPXfreeprob(cpxenv, &((*lpi)->cpxlp)) );
   freeMemory(*lpi);

   /* free environment */
   numlp--;
   if( numlp == 0 )
   {
      CHECK_ZERO( CPXcloseCPLEX(&cpxenv) );
   }

   return SCIP_OKAY;
}

RETCODE SCIPcopyLP(                     /**< copies data into LP problem object */
   LPI*             lpi,                /**< LP interface structure */
   int              ncol,               /**< number of columns */
   int              nrow,               /**< number of rows */
   OBJSEN           objsen,             /**< objective sense */
   const double*    obj,                /**< objective function vector */
   const double*    rhs,                /**< right hand side vector */
   const char*      sen,                /**< row sense vector */
   const int*       beg,                /**< start index of each column in ind- and val-array */
   const int*       cnt,                /**< number of nonzeros for each column */
   const int*       ind,                /**< row indices of constraint matrix entries */
   const double*    val,                /**< values of constraint matrix entries */
   const double*    lb,                 /**< lower bound vector */
   const double*    ub,                 /**< upper bound vector */
   const char**     cname,              /**< column names */
   const char**     rname               /**< row names */
   )
{
   int     restat;

   assert(cpxenv != NULL);
   assert(lpi != NULL);
   assert(lpi->cpxlp != NULL);

   invalidateSolution(lpi);

   CHECK_ZERO( CPXcopylpwnames(cpxenv, lpi->cpxlp, ncol, nrow, cpxObjsen(objsen), 
                  (double*)obj, (double*)rhs, (char*)sen, (int*)beg, (int*)cnt, (int*)ind, (double*)val,
                  (double*)lb, (double*)ub, NULL, (char**)cname, (char**)rname) );

   return SCIP_OKAY;
}

RETCODE SCIPaddLPRows(                  /**< adds rows to the LP */
   LPI*             lpi,                /**< LP interface structure */
   int              nrow,               /**< number of rows to be added */
   int              nnonz,              /**< number of nonzero elements to be added to the constraint matrix */
   const double*    rhs,                /**< right hand side vector of new rows */
   const char*      sen,                /**< row senses */
   const int*       beg,                /**< start index of each row in ind- and val-array */
   const int*       ind,                /**< column indices of constraint matrix entries */
   const double*    val,                /**< values of constraint matrix entries */
   const char**     name                /**< row names */
   )
{
   assert(cpxenv != NULL);
   assert(lpi != NULL);
   assert(lpi->cpxlp != NULL);

   invalidateSolution(lpi);

   CHECK_ZERO( CPXaddrows(cpxenv, lpi->cpxlp, 0, nrow, nnonz, (double*)rhs, (char*)sen,
                  (int*)beg, (int*)ind, (double*)val, NULL, (char**)name) );

   return SCIP_OKAY;
}

RETCODE SCIPdelLPRows(                  /**< deletes rows from LP */
   LPI*             lpi,                /**< LP interface structure */
   int*             dstat               /**< deletion status of rows
                                         *   input:  neg. value if row should be deleted, non-neg. value if not
                                         *   output: new position of row, -1 if row was deleted */
   )
{
   assert(cpxenv != NULL);
   assert(lpi != NULL);
   assert(lpi->cpxlp != NULL);

   invalidateSolution(lpi);

   CHECK_ZERO( CPXdelsetrows(cpxenv, lpi->cpxlp, dstat) );

   return SCIP_OKAY;   
}

RETCODE SCIPgetrowBinv(                 /**< get dense row of inverse basis matrix (A_B)^-1 */
   LPI*             lpi,                /**< LP interface structure */
   int              i,                  /**< row number */
   double*          val                 /**< vector to return coefficients */
   )
{
   assert(cpxenv != NULL);
   assert(lpi != NULL);
   assert(lpi->cpxlp != NULL);

   CHECK_ZERO( CPXbinvrow(cpxenv, lpi->cpxlp, i, val) );

   return SCIP_OKAY;
}

RETCODE SCIPgetrowBinvA(                /**< get dense row of inverse basis matrix times constraint matrix (A_B)^-1 * A */
   LPI*             lpi,                /**< LP interface structure */
   int              i,                  /**< row number */
   const double*    binv,               /**< dense row vector of row in (A_B)^-1 from prior call to SCIPgetrowBinv() */
   double*          val                 /**< vector to return coefficients */
   )
{
   assert(cpxenv != NULL);
   assert(lpi != NULL);
   assert(lpi->cpxlp != NULL);

   CHECK_ZERO( CPXbinvarow(cpxenv, lpi->cpxlp, i, val) );

   return SCIP_OKAY;
}

RETCODE SCIPgetLb(                      /**< gets lower bounds of variables */
   LPI*             lpi,                /**< LP interface structure */
   int              beg,                /**< first variable to get bound for */
   int              end,                /**< last variable to get bound for */
   double*          lb                  /**< vector to store the bounds */
   )
{
   assert(cpxenv != NULL);
   assert(lpi != NULL);
   assert(lpi->cpxlp != NULL);

   CHECK_ZERO( CPXgetlb(cpxenv, lpi->cpxlp, lb, beg, end) );

   return SCIP_OKAY;
}

RETCODE SCIPgetUb(                      /**< gets upper bounds of variables */
   LPI*             lpi,                /**< LP interface structure */
   int              beg,                /**< first variable to get bound for */
   int              end,                /**< last variable to get bound for */
   double*          ub                  /**< vector to store the bounds */
   )
{
   assert(cpxenv != NULL);
   assert(lpi != NULL);
   assert(lpi->cpxlp != NULL);

   CHECK_ZERO( CPXgetub(cpxenv, lpi->cpxlp, ub, beg, end) );

   return SCIP_OKAY;
}

RETCODE SCIPchgBd(                      /**< changes bounds of the variables in the LP */
   LPI*             lpi,                /**< LP interface structure */
   int              n,                  /**< number of bounds to be changed */
   const int*       ind,                /**< column indices */
   const char*      lu,                 /**< specifies, if 'L'ower or 'U'pper bound should be changed */
   const double*    bd                  /**< values for the new bounds */
   )
{
   assert(cpxenv != NULL);
   assert(lpi != NULL);
   assert(lpi->cpxlp != NULL);

   invalidateSolution(lpi);

   CHECK_ZERO( CPXchgbds(cpxenv, lpi->cpxlp, n, (int*)ind, (char*)lu, (double*)bd) );

   return SCIP_OKAY;
}

RETCODE SCIPchgRhs(                     /**< changes right hand sides of rows in the LP */
   LPI*             lpi,                /**< LP interface structure */
   int              n,                  /**< number of rows to change */
   const int*       ind,                /**< row indices */
   const double*    rhs                 /**< new values for right hand sides */
   )
{
   assert(cpxenv != NULL);
   assert(lpi != NULL);
   assert(lpi->cpxlp != NULL);

   invalidateSolution(lpi);

   CHECK_ZERO( CPXchgrhs(cpxenv, lpi->cpxlp, n, (int*)ind, (double*)rhs) );

   return SCIP_OKAY;
}

RETCODE SCIPchgObjsen(                  /**< changes the objective sense */
   LPI*             lpi,                /**< LP interface structure */
   OBJSEN           objsen              /**< new objective sense */
   )
{
   assert(cpxenv != NULL);
   assert(lpi != NULL);
   assert(lpi->cpxlp != NULL);

   invalidateSolution(lpi);
   
   CPXchgobjsen(cpxenv, lpi->cpxlp, cpxObjsen(objsen));

   return SCIP_OKAY;
}

RETCODE SCIPgetBind(                    /**< returns the indices of the basic columns and rows */
   LPI*             lpi,                /**< LP interface structure */
   int*             bind                /**< basic column n gives value n, basic row m gives value -1-m */
   )
{
   assert(cpxenv != NULL);
   assert(lpi != NULL);
   assert(lpi->cpxlp != NULL);

   CHECK_ZERO( CPXgetbhead(cpxenv, lpi->cpxlp, bind, NULL) );

   return SCIP_OKAY;
}

RETCODE SCIPgetLPIntpar(                /**< gets integer parameter of LP */
   LPI*             lpi,                /**< LP interface structure */
   LPPARAM          type,               /**< parameter number */
   int*             ival                /**< buffer to store the parameter value */
   )
{
   int advind;

   assert(cpxenv != NULL);
   assert(lpi != NULL);
   assert(lpi->cpxlp != NULL);
   assert(ival != NULL);

   switch( type )
   {
   case SCIP_LPPAR_FROMSCRATCH:
      if( getIntParam(lpi, CPX_PARAM_ADVIND) == CPX_ON )
	 *ival = SCIP_DISABLED;
      else
	 *ival = SCIP_ENABLED;
      break;
   case SCIP_LPPAR_LPIT1:
      *ival = CPXgetphase1cnt(cpxenv, lpi->cpxlp);
      break;
   case SCIP_LPPAR_LPIT2:
      *ival = CPXgetitcnt(cpxenv, lpi->cpxlp);
      break;
   default:
      return SCIP_LPERROR;
   }

   return SCIP_OKAY;
}

RETCODE SCIPsetLPIntpar(                /**< sets integer parameter of LP */
   LPI*             lpi,                /**< LP interface structure */
   LPPARAM          type,               /**< parameter number */
   int              ival                /**< parameter value */
   )
{
   assert(cpxenv != NULL);
   assert(lpi != NULL);
   assert(lpi->cpxlp != NULL);

   switch( type )
   {
   case SCIP_LPPAR_FROMSCRATCH:
      assert(ival == SCIP_ENABLED || ival == SCIP_DISABLED);
      setIntParam(lpi, CPX_PARAM_ADVIND, (ival == SCIP_DISABLED) ? CPX_ON : CPX_OFF);
      break;
   case SCIP_LPPAR_LPITLIM:
      setIntParam(lpi, CPX_PARAM_ITLIM, ival);
      break;
   case SCIP_LPPAR_FASTMIP:
      assert(ival == SCIP_ENABLED || ival == SCIP_DISABLED);
      setIntParam(lpi, CPX_PARAM_FASTMIP, (ival == SCIP_ENABLED) ? CPX_ON : CPX_OFF);
      break;
   case SCIP_LPPAR_PRICING:
      switch( (PRICING)ival )
      {
      case SCIP_PRICING_FULL:
	 setIntParam(lpi, CPX_PARAM_DPRIIND, CPX_DPRIIND_FULL);
         break;
      case SCIP_PRICING_STEEP:
	 setIntParam(lpi, CPX_PARAM_DPRIIND, CPX_DPRIIND_STEEP);
	 break;
      case SCIP_PRICING_STEEPQSTART:
	 setIntParam(lpi, CPX_PARAM_DPRIIND, CPX_DPRIIND_STEEPQSTART);
	 break;
      default:
         return SCIP_LPERROR;
      }
      break;
   case SCIP_LPPAR_LPINFO:
      assert(ival == SCIP_ENABLED || ival == SCIP_DISABLED);
      if( ival == SCIP_ENABLED )
      {
	 setIntParam(lpi, CPX_PARAM_SIMDISPLAY, CPX_ON);
	 setIntParam(lpi, CPX_PARAM_SCRIND, CPX_ON);
      }
      else 
      {
	 setIntParam(lpi, CPX_PARAM_SIMDISPLAY, CPX_OFF);
	 setIntParam(lpi, CPX_PARAM_SCRIND, CPX_OFF);
      }
      break;
   default:
      return SCIP_LPERROR;
   }

   return SCIP_OKAY;
}

RETCODE SCIPgetLPDblpar(                /**< gets double parameter of LP */
   LPI*             lpi,                /**< LP interface structure */
   LPPARAM          type,               /**< parameter number */
   double*          dval                /**< buffer to store the parameter value */
   )
{
   assert(cpxenv != NULL);
   assert(lpi != NULL);
   assert(lpi->cpxlp != NULL);
   assert(dval != NULL);

   switch( type )
   {
   case SCIP_LPPAR_FEASTOL:
      *dval = getDblParam(lpi, CPX_PARAM_EPRHS);
      break;
   default:
      return SCIP_LPERROR;
      break;
   }
   
   return SCIP_OKAY;
}

RETCODE SCIPsetLPDblpar(                /**< sets double parameter of LP */
   LPI*             lpi,                /**< LP interface structure */
   LPPARAM          type,               /**< parameter number */
   double           dval                /**< parameter value */
   )
{
   assert(cpxenv != NULL);
   assert(lpi != NULL);
   assert(lpi->cpxlp != NULL);

   switch( type )
   {
   case SCIP_LPPAR_FEASTOL:
      setDblParam(lpi, CPX_PARAM_EPRHS, dval);
      break;
   case SCIP_LPPAR_LOBJLIM:
      setDblParam(lpi, CPX_PARAM_OBJLLIM, dval);
      break;
   case SCIP_LPPAR_UOBJLIM:
      setDblParam(lpi, CPX_PARAM_OBJULIM, dval);
      break;
   case SCIP_LPPAR_LPTILIM:
      setDblParam(lpi, CPX_PARAM_TILIM, dval);
      break;
   default:
      return SCIP_LPERROR;
   }

   return SCIP_OKAY;
}

RETCODE SCIPgetSol(                     /**< gets primal and dual solution vectors */
   LPI*             lpi,                /**< LP interface structure */
   double*          objval,             /**< stores the objective value */
   double*          psol,               /**< primal solution vector */
   double*          pi,                 /**< dual solution vector */
   double*          slck,               /**< slack vector */
   double*          redcost             /**< reduced cost vector */
   )
{
   int dummy;

   assert(cpxenv != NULL);
   assert(lpi != NULL);
   assert(lpi->cpxlp != NULL);

   CHECK_ZERO( CPXsolution(cpxenv, lpi->cpxlp, &dummy, objval, psol, pi, slck, redcost) );
   assert(dummy == lpi->solstat);

   return SCIP_OKAY;
}

RETCODE SCIPstrongbranch(               /**< performs strong branching iterations on all candidates */
   LPI*             lpi,                /**< LP interface structure */
   const double*    psol,               /**< primal LP solution vector */
   int              ncand,              /**< size of candidate list */
   const int*       cand,               /**< candidate list */
   int              itlim,              /**< iteration limit for strong branchings */
   double*          down,               /**< stores dual bound after branching candidate down */
   double*          up                  /**< stores dual bound after branching candidate up */
   )
{
   assert(cpxenv != NULL);
   assert(lpi != NULL);
   assert(lpi->cpxlp != NULL);

   CHECK_ZERO( CPXstrongbranch(cpxenv, lpi->cpxlp, (int*)cand, ncand, down, up, itlim) );

   return SCIP_OKAY;
}

RETCODE SCIPoptLPPrimal(                /**< calls primal simplex to solve the LP */
   LPI*             lpi                 /**< LP interface structure */
   )
{
   assert(cpxenv != NULL);
   assert(lpi != NULL);
   assert(lpi->cpxlp != NULL);

   invalidateSolution(lpi);

   setParameterValues(&(lpi->cpxparam));

   lpi->optalgo = PRIMALSIMPLEX;
   lpi->retval = CPXprimopt( cpxenv, lpi->cpxlp );
   switch( lpi->retval  )
   {
   case 0:
   case CPXERR_PRESLV_INForUNBD:
   case CPXERR_PRESLV_INF:
   case CPXERR_PRESLV_UNBD:
      break;
   case CPXERR_NO_MEMORY:
      return SCIP_NOMEMORY;
   default:
      return SCIP_LPERROR;
   }

   lpi->solstat = CPXgetstat(cpxenv, lpi->cpxlp);
   
   return SCIP_OKAY;
}

RETCODE SCIPoptLPDual(                  /**< calls dual simplex to solve the LP */
   LPI*             lpi                 /**< LP interface structure */
   )
{
   assert(cpxenv != NULL);
   assert(lpi != NULL);
   assert(lpi->cpxlp != NULL);

   invalidateSolution(lpi);

   setParameterValues(&(lpi->cpxparam));

   lpi->optalgo = DUALSIMPLEX;
   lpi->retval = CPXdualopt( cpxenv, lpi->cpxlp );
   switch( lpi->retval  )
   {
   case 0:
   case CPXERR_PRESLV_INForUNBD:
   case CPXERR_PRESLV_INF:
   case CPXERR_PRESLV_UNBD:
      break;
   case CPXERR_NO_MEMORY:
      return SCIP_NOMEMORY;
   default:
      return SCIP_LPERROR;
   }

   lpi->solstat = CPXgetstat(cpxenv, lpi->cpxlp);
   
   return SCIP_OKAY;
}

Bool SCIPisPrimalUnbounded(             /**< returns TRUE iff LP is primal unbounded */
   LPI*             lpi                 /**< LP interface structure */
   )
{
   assert(cpxenv != NULL);
   assert(lpi != NULL);
   assert(lpi->cpxlp != NULL);

   switch( lpi->optalgo )
   {
   case PRIMALSIMPLEX:
      return( isUnboundedSolution(lpi) || lpi->solstat == CPX_UNBOUNDED );
   case DUALSIMPLEX:
      /* primal unbounded means dual infeasible */
      return( isInfeasibleSolution(lpi) || lpi->solstat == CPX_INFEASIBLE );
   default:
      errorMessage("LP not optimized");
      return FALSE;
   }
}

Bool SCIPisPrimalInfeasible(            /**< returns TRUE iff LP is primal infeasible */
   LPI*             lpi                 /**< LP interface structure */
   )
{
   assert(cpxenv != NULL);
   assert(lpi != NULL);
   assert(lpi->cpxlp != NULL);

   switch( lpi->optalgo )
   {
   case PRIMALSIMPLEX:
      return( isInfeasibleSolution(lpi) || lpi->solstat == CPX_INFEASIBLE );
   case DUALSIMPLEX:
      /* primal infeasible means dual unbounded */
      return( isUnboundedSolution(lpi) || lpi->solstat == CPX_UNBOUNDED );
   default:
      errorMessage("LP not optimized");
      return FALSE;
   }
}

Bool SCIPisOptimal(                     /**< returns TRUE iff LP was solved to optimality */
   LPI*             lpi                 /**< LP interface structure */
   )
{
   assert(cpxenv != NULL);
   assert(lpi != NULL);
   assert(lpi->cpxlp != NULL);

   return( isOptimalSolution(lpi) && lpi->solstat == CPX_OPTIMAL );
}

Bool SCIPisDualValid(                   /**< returns TRUE iff actual LP solution is dual valid */
   LPI*             lpi                 /**< LP interface structure */
   )
{
   assert(cpxenv != NULL);
   assert(lpi != NULL);
   assert(lpi->cpxlp != NULL);

   return( isValidSolution(lpi)
      && ( lpi->solstat == CPX_OPTIMAL
	 || lpi->solstat == CPX_OBJ_LIM
         || lpi->solstat == CPX_IT_LIM_FEAS
	 || lpi->solstat == CPX_TIME_LIM_FEAS
         || lpi->solstat == CPX_NUM_BEST_FEAS
	 || lpi->solstat == CPX_ABORT_FEAS ) );
}

Bool SCIPisStable(                      /**< returns TRUE iff actual LP basis is stable */
   LPI*             lpi                 /**< LP interface structure */
   )
{
   assert(cpxenv != NULL);
   assert(lpi != NULL);
   assert(lpi->cpxlp != NULL);

   return( isValidSolution(lpi)
      && lpi->solstat != CPX_NUM_BEST_FEAS
      && lpi->solstat != CPX_NUM_BEST_INFEAS
      && lpi->solstat != CPX_OPTIMAL_INFEAS );
}

Bool SCIPisLPError(                     /**< returns TRUE iff an error occured while solving the LP */
   LPI*             lpi                 /**< LP interface structure */
   )
{
   assert(cpxenv != NULL);
   assert(lpi != NULL);
   assert(lpi->cpxlp != NULL);

   return !isValidSolution(lpi);
}

Bool SCIPisObjlimExc(                   /**< returns TRUE iff the objective limit was reached */
   LPI*             lpi                 /**< LP interface structure */
   )
{
   assert(cpxenv != NULL);
   assert(lpi != NULL);
   assert(lpi->cpxlp != NULL);

   return( isUnboundedSolution(lpi)
      || lpi->solstat == CPX_UNBOUNDED
      || lpi->solstat == CPX_OBJ_LIM );
}

Bool SCIPisIterlimExc(                  /**< returns TRUE iff the iteration limit was reached */
   LPI*             lpi                 /**< LP interface structure */
   )
{
   assert(cpxenv != NULL);
   assert(lpi != NULL);
   assert(lpi->cpxlp != NULL);

   return( lpi->solstat == CPX_IT_LIM_FEAS
      || lpi->solstat == CPX_IT_LIM_INFEAS );
}

Bool SCIPisTimelimExc(                  /**< returns TRUE iff the time limit was reached */
   LPI*             lpi                 /**< LP interface structure */
   )
{
   assert(cpxenv != NULL);
   assert(lpi != NULL);
   assert(lpi->cpxlp != NULL);

   return( lpi->solstat == CPX_TIME_LIM_FEAS
      || lpi->solstat == CPX_TIME_LIM_INFEAS );
}


static 
int colpacketNum(                       /**< returns the number of packets needed to store column packet information */
   int              ncol                /**< number of columns to store */
   )
{
   return (ncol+COLS_PER_PACKET-1)/COLS_PER_PACKET;
}

static 
int rowpacketNum(                       /**< returns the number of packets needed to store row packet information */
   int              nrow                /**< number of rows to store */
   )
{
   return (nrow+ROWS_PER_PACKET-1)/ROWS_PER_PACKET;
}

static
void packLPState(                       /**< store row and column basis status in a packed LP state object */
   LPSTATE*        lpstate,             /**< pointer to LP state data */
   const int*      cstat,               /**< basis status of columns in unpacked format */
   const int*      rstat                /**< basis status of rows in unpacked format */
   )
{
   assert(lpstate != NULL);
   assert(lpstate->packcstat != NULL);
   assert(lpstate->packrstat != NULL);

   SCIPencodeDualBit(cstat, lpstate->packcstat, lpstate->ncol);
   SCIPencodeSingleBit(rstat, lpstate->packrstat, lpstate->nrow);
}

static
void unpackLPState(                     /**< unpack row and column basis status from a packed LP state object */
   const LPSTATE*  lpstate,             /**< pointer to LP state data */
   int*            cstat,               /**< buffer for storing basis status of columns in unpacked format */
   int*            rstat                /**< buffer for storing basis status of rows in unpacked format */
   )
{
   assert(lpstate != NULL);
   assert(lpstate->packcstat != NULL);
   assert(lpstate->packrstat != NULL);

   SCIPdecodeDualBit(lpstate->packcstat, cstat, lpstate->ncol);
   SCIPdecodeSingleBit(lpstate->packrstat, rstat, lpstate->nrow);
}

static
LPSTATE* createLPState(
   MEM*             mem,                /**< block memory buffers */
   int              ncol,               /**< number of columns to store */
   int              nrow                /**< number of rows to store */
   )
{
   LPSTATE* lpstate;

   assert(mem != NULL);
   assert(ncol >= 0);
   assert(nrow >= 0);

   ALLOC_NULL( allocBlockMemory(mem->statemem, lpstate) );
   ALLOC_NULL( allocBlockMemoryArray(mem->statemem, lpstate->packcstat, colpacketNum(ncol)) );
   ALLOC_NULL( allocBlockMemoryArray(mem->statemem, lpstate->packrstat, rowpacketNum(nrow)) );
   lpstate->dnorm = NULL;

   return lpstate;
}

static
void freeLPState(
   MEM*             mem,                /**< block memory buffers */
   LPSTATE**        lpstate             /**< pointer to LP state information (like basis information) */
   )
{
   assert(mem != NULL);
   assert(lpstate != NULL);
   assert(*lpstate != NULL);
   assert((*lpstate)->numuses == 0);

   freeBlockMemoryArray(mem->statemem, (*lpstate)->packcstat, colpacketNum((*lpstate)->ncol));
   freeBlockMemoryArray(mem->statemem, (*lpstate)->packrstat, rowpacketNum((*lpstate)->nrow));
   freeBlockMemoryArrayNull(mem->statemem, (*lpstate)->dnorm, (*lpstate)->ncol);
   freeBlockMemory(mem->statemem, *lpstate);
}

RETCODE SCIPsaveLPState(                /**< stores LP state (like basis information) */
   MEM*             mem,                /**< block memory buffers */
   LPI*             lpi,                /**< LP interface structure */
   LPSTATE**        lpstate             /**< pointer to LP state information (like basis information) */
   )
{
   int  ncol;
   int  nrow;
   int* cstat;
   int* rstat;

   assert(mem != NULL);
   assert(cpxenv != NULL);
   assert(lpi != NULL);
   assert(lpi->cpxlp != NULL);
   assert(lpstate != NULL);

   ncol = CPXgetnumcols(cpxenv, lpi->cpxlp);
   nrow = CPXgetnumrows(cpxenv, lpi->cpxlp);
   assert(0 <= ncol && ncol < SCIP_MAXNCOL);
   assert(0 <= nrow && nrow < SCIP_MAXNROW);
   
   /* allocate lpstate data */
   ALLOC_OKAY( *lpstate = createLPState(mem, ncol, nrow) );

   /* allocate temporary buffer for storing uncompressed basis information */
   ALLOC_OKAY( allocBlockMemoryArray(mem->tempmem, cstat, ncol) );
   ALLOC_OKAY( allocBlockMemoryArray(mem->tempmem, rstat, nrow) );

   if( getIntParam(lpi, CPX_PARAM_DPRIIND) == CPX_DPRIIND_STEEP )
   {
      ALLOC_OKAY( allocBlockMemoryArray(mem->statemem, (*lpstate)->dnorm, ncol) );
      CHECK_ZERO( CPXgetbasednorms(cpxenv, lpi->cpxlp, cstat, rstat, (*lpstate)->dnorm) );
   }
   else
   {
      (*lpstate)->dnorm = NULL;
      CHECK_ZERO( CPXgetbase(cpxenv, lpi->cpxlp, cstat, rstat) );
   }

   /* fill LP state data */
   (*lpstate)->ncol = ncol;
   (*lpstate)->nrow = nrow;
   packLPState(*lpstate, cstat, rstat);

   /* free temporary memory */
   freeBlockMemoryArray(mem->tempmem, cstat, ncol);
   freeBlockMemoryArray(mem->tempmem, rstat, nrow);

   return SCIP_OKAY;
}

RETCODE SCIPloadLPState(                /**< loads LP state (like basis information) into solver */
   MEM*             mem,                /**< block memory buffers */
   LPI*             lpi,                /**< LP interface structure */
   LPSTATE*         lpstate             /**< LP state information (like basis information) */
   )
{
   int  pricing;
   int* cstat;
   int* rstat;

   assert(mem != NULL);
   assert(cpxenv != NULL);
   assert(lpi != NULL);
   assert(lpi->cpxlp != NULL);
   assert(lpstate != NULL);
   assert(lpstate->ncol == CPXgetnumcols(cpxenv, lpi->cpxlp));
   assert(lpstate->nrow == CPXgetnumrows(cpxenv, lpi->cpxlp));

   /* allocate temporary buffer for storing uncompressed basis information */
   ALLOC_OKAY( allocBlockMemoryArray(mem->tempmem, cstat, lpstate->ncol) );
   ALLOC_OKAY( allocBlockMemoryArray(mem->tempmem, rstat, lpstate->nrow) );

   unpackLPState(lpstate, cstat, rstat);
   if( getIntParam(lpi, CPX_PARAM_DPRIIND) == CPX_DPRIIND_STEEP && lpstate->dnorm != NULL )
   {
      CHECK_ZERO( CPXcopybasednorms(cpxenv, lpi->cpxlp, cstat, rstat, lpstate->dnorm) );
   }
   else
   {
      CHECK_ZERO( CPXcopybase(cpxenv, lpi->cpxlp, cstat, rstat) );
   }

   /* free temporary memory */
   freeBlockMemoryArray(mem->tempmem, cstat, lpstate->ncol);
   freeBlockMemoryArray(mem->tempmem, rstat, lpstate->nrow);

   return SCIP_OKAY;
}

void SCIPcaptureLPState(                /**< increases usage counter of LP state */
   LPSTATE*         lpstate             /**< LP state information (like basis information) */
   )
{
   assert(lpstate != NULL);
   assert(lpstate->numuses >= 0);

   lpstate->numuses++;
}

void SCIPreleaseLPState(                /**< decreases usage counter of LP state, and frees memory if necessary */
   MEM*             mem,                /**< block memory buffers */
   LPSTATE**        lpstate             /**< LP state information (like basis information) */
   )
{
   assert(mem != NULL);
   assert(lpstate != NULL);
   assert(*lpstate != NULL);
   assert((*lpstate)->numuses >= 1);

   (*lpstate)->numuses--;
   if( (*lpstate)->numuses == 0 )
      freeLPState(mem, lpstate);
}

RETCODE SCIPwriteLPState(               /**< writes LP state (like basis information) to a file */
   LPI*             lpi,                /**< LP interface structure */
   const char*      fname               /**< file name */
   )
{
   assert(cpxenv != NULL);
   assert(lpi != NULL);
   assert(lpi->cpxlp != NULL);

   CHECK_ZERO( CPXmbasewrite(cpxenv, lpi->cpxlp, (char*)fname) );

   return SCIP_OKAY;
}

RETCODE SCIPwriteLP(                    /**< writes LP to a file */
   LPI*             lpi,                /**< LP interface structure */
   const char*      fname               /**< file name */
   )
{
   assert(cpxenv != NULL);
   assert(lpi != NULL);
   assert(lpi->cpxlp != NULL);

   CHECK_ZERO( CPXlpwrite(cpxenv, lpi->cpxlp, (char*)fname) );

   return SCIP_OKAY;
}
