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

/**@file   lpi_spx.cpp
 * @ingroup LPIS
 * @brief  LP interface for SoPlex version 1.4 and higher
 * @author Tobias Achterberg
 * @author Timo Berthold
 * @author Ambros Gleixner
 * @author Marc Pfetsch
 *
 * This is an implementation of SCIP's LP interface for SoPlex. While the ratio test is fixed to SoPlex's standard,
 * different pricing methods can be chosen and an autopricing strategy (start with devex and switch to steepest edge
 * after too many iterations) is implemented directly. Scaler and simplifier may be applied if solving from scratch.
 *
 * For debugging purposes, the SoPlex results can be double checked with CPLEX if WITH_LPSCHECK is defined. This may
 * yield false positives, since the LP is dumped to a file for transfering it to CPLEX, hence, precision may be lost.
 */

/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#define AUTOPRICING_ITERSWITCH          10000/**< start with devex and switch to steepest edge after this many iterations */
#define STRONGBRANCH_RESTOREBASIS            /**< if defined then in SCIPlpiStrongbranch() we restore the basis after the
                                              *   down branch and after the up branch; if false only after the end of a
                                              *   strong branching phase, which however seems to mostly increase strong
                                              *   branching time and iterations */

/* in this case the SoPlex results are double checked using CPLEX */
#ifdef WITH_LPSCHECK
#include <cplex.h>

#define CHECK_SPXSOLVE                  true /**< shall the SoPlex results in spxSolve() be double checked using CPLEX? */
#define CHECK_SPXSTRONGBRANCH           true /**< shall the SoPlex results in SCIPlpStrongbranch() be double checked using CPLEX? */
#define CHECK_START                     0    /**< skip first CHECK_START number of checks */
#define EXIT_AT_WRONG_RESULT            false/**< shall program be exited if CPLEX returns different result than SoPlex? */
#define EXIT_AT_CPXERROR                false/**< shall program be exited if CPLEX returns an error? */

#define CPX_CALL(x)                     do                                                                                  \
                                        {                                                                                   \
                                           int _cpxstat_;               \
                                           if( (_cpxstat_ = (x)) != 0 )                                                     \
                                           {                                                                                \
                                              SCIPmessagePrintWarning(m_messagehdlr, "CPLEX error <%d>; SoPlex result unchecked\n", _cpxstat_); \
                                              if( EXIT_AT_CPXERROR )                                                        \
                                              {                                                                             \
                                                 exit(1);                                                                   \
                                              }                                                                             \
                                              else                                                                          \
                                              {                                                                             \
                                                 goto ENDCHECK;                                                             \
                                              }                                                                             \
                                           }                                                                                \
                                        }                                                                                   \
                                        while( false )
#endif

/* remember the original value of the SCIP_DEBUG define and undefine it */
#ifdef SCIP_DEBUG
#define ___DEBUG
#undef SCIP_DEBUG
#endif

/* include SoPlex solver */
#include "spxsolver.h"

/* define subversion for versions <= 1.5.0.1 */
#ifndef SOPLEX_SUBVERSION
#define SOPLEX_SUBVERSION 0
#endif

/* check version */
#if (SOPLEX_VERSION < 133)
#error "This interface is not compatible with SoPlex versions prior to 1.4"
#endif

/* get githash of SoPlex */
#if (SOPLEX_VERSION >= 160)
#include "spxgithash.h"
#endif

/* include SoPlex components */
#include "slufactor.h"
#include "spxsteeppr.h"
#if ((SOPLEX_VERSION == 160 && SOPLEX_SUBVERSION >= 6) || SOPLEX_VERSION > 160)
#include "spxsteepexpr.h"
#endif
#include "spxparmultpr.h"
#include "spxdevexpr.h"
#include "spxfastrt.h"
#include "spxmainsm.h"
#include "spxequilisc.h"

#ifdef WITH_BOUNDFLIPPING
#include "spxboundflippingrt.h"
#endif

/* reset the SCIP_DEBUG define to its original SCIP value */
#undef SCIP_DEBUG
#ifdef ___DEBUG
#define SCIP_DEBUG
#undef ___DEBUG
#endif

#define SOPLEX_VERBLEVEL                5    /**< verbosity level for LPINFO */

#include "scip/pub_message.h"

/********************************************************************/
/*----------------------------- C++ --------------------------------*/
/********************************************************************/

/* in C++ we have to use "0" instead of "(void*)0" */
#undef NULL
#define NULL 0

#include <cassert>
using namespace soplex;


/** Macro for a single SoPlex call for which exceptions have to be catched - return an LP error. We
 *  make no distinction between different exception types, e.g., between memory allocation and other
 *  exceptions.
 */
#ifndef NDEBUG
#define SOPLEX_TRY(messagehdlr, x)  do                                  \
   {                                                                    \
      try                                                               \
      {                                                                 \
         (x);                                                           \
      }                                                                 \
      catch(SPxException E)                                             \
      {                                                                 \
         std::string s = E.what();                                      \
         SCIPmessagePrintWarning((messagehdlr), "SoPlex threw an exception: %s\n", s.c_str()); \
         return SCIP_LPERROR;                                           \
      }                                                                 \
   }                                                                    \
   while( FALSE )

#else
#define SOPLEX_TRY(messagehdlr, x)  do                                  \
   {                                                                    \
      try                                                               \
      {                                                                 \
         (x);                                                           \
      }                                                                 \
      catch(SPxException E)                                             \
      {                                                                 \
         return SCIP_LPERROR;                                           \
      }                                                                 \
   }                                                                    \
   while( FALSE )
#endif

#define SOPLEX_TRYLPI(x) SOPLEX_TRY(lpi->messagehdlr, x)
#define SOPLEX_TRYLPIPTR(x) SOPLEX_TRY((*lpi)->messagehdlr, x)

/* Macro for a single SoPlex call for which exceptions have to be catched - abort if they
 * arise. SCIP_ABORT() is not accessible here.
 */
#define SOPLEX_TRY_ABORT(x)  do                                         \
   {                                                                    \
      try                                                               \
      {                                                                 \
         (x);                                                           \
      }                                                                 \
      catch(SPxException E)                                             \
      {                                                                 \
         std::string s = E.what();                                      \
         SCIPerrorMessage("SoPlex threw an exception: %s\n", s.c_str()); \
         abort();                                                       \
      }                                                                 \
   }                                                                    \
   while( FALSE )



/** SCIP's SoPlex class */
class SPxSCIP : public SPxSolver
{
   SPxLP::SPxSense       m_sense;            /**< optimization sense */
   SLUFactor             m_slu;              /**< sparse LU factorization */
   SPxSteepPR            m_price_steep;      /**< steepest edge pricer */
#if ((SOPLEX_VERSION == 160 && SOPLEX_SUBVERSION >= 6) || SOPLEX_VERSION > 160)
   SPxSteepExPR          m_price_steep_ex;   /**< steepest edge with exact weight initialization */
#else
   SPxSteepPR            m_price_steep_ex;   /**< fallback to quick start pricer */
#endif
   SPxParMultPR          m_price_parmult;    /**< partial multiple pricer */
   SPxDevexPR            m_price_devex;      /**< devex pricer */
#ifdef WITH_BOUNDFLIPPING
   SPxBoundFlippingRT    m_ratio;            /**< Long step dual ratio tester */
#else
   SPxFastRT             m_ratio;            /**< Harris fast ratio tester */
#endif
   char*                 m_probname;         /**< problem name */
   bool                  m_fromscratch;      /**< use old basis indicator */
   bool                  m_scaling;          /**< use lp scaling */
   bool                  m_presolving;       /**< use lp presolving */
   Real                  m_lpifeastol;       /**< feastol set by SCIPlpiSetRealpar() */
   Real                  m_lpiopttol;        /**< opttol set by SCIPlpiSetRealpar() */
   Real                  m_objLoLimit;       /**< lower objective limit */
   Real                  m_objUpLimit;       /**< upper objective limit */
   Status                m_stat;             /**< solving status */
   bool                  m_lpinfo;           /**< storing whether output is turned on */
   bool                  m_autopricing;      /**< is automatic pricing selected? */
   int                   m_itlim;            /**< iteration limit (-1 for unbounded) */
   int                   m_itused;           /**< number of iterations spent in phase one of auto pricing */
   SPxSolver::VarStatus* m_rowstat;          /**< basis status of rows before starting strong branching (if available, 0 otherwise) */
   SPxSolver::VarStatus* m_colstat;          /**< basis status of columns before starting strong branching (if available, 0 otherwise) */
   NameSet*              m_rownames;         /**< row names */
   NameSet*              m_colnames;         /**< column names */

#ifdef WITH_LPSCHECK
   int                   m_checknum;
   bool                  m_doublecheck;
   CPXENVptr             m_cpxenv;           /**< CPLEX memory environment */
   CPXLPptr              m_cpxlp;            /**< CPLEX lp structure */
#endif
   SCIP_MESSAGEHDLR*     m_messagehdlr;      /**< messagehdlr handler to printing messages, or NULL */

public:
   SPxSCIP(
      SCIP_MESSAGEHDLR*  messagehdlr = NULL, /**< message handler */
      const char*        probname = NULL     /**< name of problem */
      )
      : SPxSolver(LEAVE, COLUMN),
        m_probname(0),
        m_fromscratch(false),
        m_scaling(true),
        m_presolving(true),
        m_objLoLimit(-soplex::infinity),
        m_objUpLimit(soplex::infinity),
        m_stat(NO_PROBLEM),
        m_lpinfo(false),
        m_autopricing(true),
        m_itlim(-1),
        m_itused(0),
        m_rowstat(NULL),
        m_colstat(NULL),
        m_rownames(0),
        m_colnames(0),
        m_messagehdlr(messagehdlr)
   {
      m_sense = sense();
      setSense(SPxLP::MINIMIZE);
      setSolver(&m_slu);
      setTester(&m_ratio);
      setPricer(&m_price_steep);
      /* no starter */

      if ( probname != NULL )
         SOPLEX_TRY_ABORT( setProbname(probname) );

#if ((SOPLEX_VERSION == 160 && SOPLEX_SUBVERSION >= 5) || SOPLEX_VERSION > 160)
      m_lpifeastol = SPxSolver::feastol();
      m_lpiopttol = SPxSolver::opttol();
#else
      m_lpifeastol = SPxSolver::delta();
      m_lpiopttol = SPxSolver::delta();
#endif

#ifdef WITH_LPSCHECK
      int cpxstat;
      m_cpxenv = CPXopenCPLEX(&cpxstat);
      assert(m_cpxenv != NULL);
      m_cpxlp = CPXcreateprob(m_cpxenv, &cpxstat, probname != NULL ? probname : "spxcheck");
      (void) CPXsetintparam(m_cpxenv, CPX_PARAM_SCRIND, 0);
      m_checknum = 0;
#endif
   }

   virtual ~SPxSCIP()
   {
      if( m_probname != NULL )
         spx_free(m_probname);  /*lint !e1551*/

      freePreStrongbranchingBasis();

#ifdef WITH_LPSCHECK
      (void) CPXfreeprob(m_cpxenv, &m_cpxlp);
      (void) CPXcloseCPLEX(&m_cpxenv);
#endif
   }

   /**< return feastol set by SCIPlpiSetRealpar(), which might be tighter than what SoPlex accepted */
   Real feastol()
   {
      return m_lpifeastol;
   }

   /**< set feastol and store value in case SoPlex only accepts a larger tolerance */
   void setFeastol(
      const Real d
      )
   {
      m_lpifeastol = d;

#if ((SOPLEX_VERSION == 160 && SOPLEX_SUBVERSION >= 5) || SOPLEX_VERSION > 160)
      SPxSolver::setFeastol(d);
#else
      SPxSolver::setDelta(d);
#endif
   }

   /**< return opttol set by SCIPlpiSetRealpar(), which might be tighter than what SoPlex accepted */
   Real opttol()
   {
      return m_lpiopttol;
   }

   /**< set opttol and store value in case SoPlex only accepts a larger tolerance */
   void setOpttol(
      const Real d
      )
   {
      m_lpiopttol = d;

#if ((SOPLEX_VERSION == 160 && SOPLEX_SUBVERSION >= 5) || SOPLEX_VERSION > 160)
      SPxSolver::setOpttol(d);
#else
      SPxSolver::setDelta(d);
#endif
   }

   bool isPerturbed()
   {
      /* the epsilon is 1e-16; we add a factor of ten to account for numerics */
      return (shift() >= 10.0 * epsilon());
   }

   /** set iteration limit (-1 = unbounded) */
   void setIterationLimit(
      const int          itlim
      )
   {
      m_itlim = itlim;
   }

   void setAutoPricer()
   {
      setPricer(&m_price_devex);
      m_autopricing = true;
   }

   void setFullPricer()
   {
      setPricer(&m_price_steep);
      m_autopricing = false;
   }

   void setSteepPricer()
   {
      setPricer(&m_price_steep_ex);
      m_autopricing = false;
   }

   void setSteepQStartPricer()
   {
      setPricer(&m_price_steep);
      m_autopricing = false;
   }

   void setParmultPricer()
   {
      setPricer(&m_price_parmult);
      m_autopricing = false;
   }

   void setDevexPricer()
   {
      setPricer(&m_price_devex);
      m_autopricing = false;
   }

   /** get iteration limit (-1 = unbounded) */
   int getIterationLimit()
   {
      return m_itlim;
   }

   bool getFromScratch() const
   {
      return m_fromscratch;
   }

   void setFromScratch(bool fs)
   {
      m_fromscratch = fs;
   }

   bool getScaling() const
   {
      return m_scaling;
   }

   void setScaling(bool s)
   {
      m_scaling = s;
   }

   bool getPresolving() const
   {
      return m_presolving;
   }

   void setPresolving(bool p)
   {
      m_presolving = p;
   }

   bool getLpInfo() const
   {
      return m_lpinfo;
   }

   void setLpInfo(bool li)
   {
      m_lpinfo = li;
   }

   SPxLP::SPxSense getSense() const
   {
      assert(m_sense == sense());

      return m_sense;
   }

   void setSense(const SPxLP::SPxSense sen)
   {
      assert(m_sense == sense());

      if( m_sense != sen )
      {
         m_sense = sen;
         changeSense(sen);

         /* if objective limit was set for the new sense previously, we have to apply it now */
         if( m_sense == SPxLP::MINIMIZE && getObjUpLimit() < soplex::infinity )
         {
            SCIPdebugMessage("setting termination value to <%g>\n", getObjUpLimit());
            SPxSolver::setTerminationValue(getObjUpLimit());
         }
         else if( m_sense == SPxLP::MAXIMIZE && getObjLoLimit() > -soplex::infinity )
         {
            SCIPdebugMessage("setting termination value to <%g>\n", getObjLoLimit());
            SPxSolver::setTerminationValue(getObjLoLimit());
         }
      }
   }

   void setProbname(const char* probname)
   {
      int len;

      assert(probname != NULL);
      if( m_probname != NULL )
         spx_free(m_probname);
      len = (int)strlen(probname);
      spx_alloc(m_probname, len + 1);
      strncpy(m_probname, probname, len);
      m_probname[len] = '\0';
   }

   Real getObjLoLimit() const
   {
      return m_objLoLimit;
   }

   void setObjLoLimit(Real limit)
   {
      if( getSense() == SPxLP::MAXIMIZE )
      {
         SCIPdebugMessage("setting termination value from <%g> to <%g>\n", m_objLoLimit, limit);
         SPxSolver::setTerminationValue(limit);
      }
      m_objLoLimit = limit;
   }

   Real getObjUpLimit() const
   {
      return m_objUpLimit;
   }

   void setObjUpLimit(Real limit)
   {
      if( getSense() == SPxLP::MINIMIZE )
      {
         SCIPdebugMessage("setting termination value from <%g> to <%g>\n", m_objUpLimit, limit);
         SPxSolver::setTerminationValue(limit);
      }
      m_objUpLimit = limit;
   }

   void setRep(SPxSolver::Representation p_rep)
   {
      if( p_rep != rep() )
      {
         SCIPdebugMessage("switching to %s representation of the basis\n", p_rep == SPxSolver::ROW ? "row" : "column");
         SPxSolver::setRep(p_rep);
      }
   }

#ifdef WITH_LPSCHECK
   bool getDoubleCheck()
   {
      m_checknum++;
      return m_doublecheck && m_checknum + 1 >= CHECK_START;
   }

   void setDoubleCheck(bool dc)
   {
      m_doublecheck = dc;
   }

   const char* spxStatusString(const SPxSolver::Status stat)
   {
      switch( stat )
      {
      case SPxSolver::ABORT_TIME:
         return "ABORT_TIME";
      case SPxSolver::ABORT_ITER:
         return "ABORT_ITER";
      case SPxSolver::ABORT_VALUE:
         return "ABORT_VALUE";
      case SPxSolver::SINGULAR:
         return "SINGULAR";
      case SPxSolver::REGULAR:
         return "REGULAR";
      case SPxSolver::UNKNOWN:
         return "UNKNOWN";
      case SPxSolver::OPTIMAL:
         return "OPTIMAL";
      case SPxSolver::UNBOUNDED:
         return "UNBOUNDED";
      case SPxSolver::INFEASIBLE:
         return "INFEASIBLE";
      default:
         return "UNKNOWN";
      }  /*lint !e788*/

      return "UNKNOWN";
   }

   const char* cpxStatusString(const int stat)
   {
      switch( stat )
      {
      case CPX_STAT_ABORT_TIME_LIM:
         return "ABORT_TIME";
      case CPX_STAT_ABORT_IT_LIM:
         return "ABORT_ITER";
      case CPX_STAT_ABORT_OBJ_LIM:
         return "ABORT_VALUE";
      case CPX_STAT_OPTIMAL:
         return "OPTIMAL";
      case CPX_STAT_OPTIMAL_INFEAS:
         return "CPX_STAT_OPTIMAL_INFEAS: OPT SOL INFEASIBLE AFTER UNSCALING";
      case CPX_STAT_UNBOUNDED:
         return "UNBOUNDED";
      case CPX_STAT_INFEASIBLE:
         return "INFEASIBLE";
      case CPX_STAT_INForUNBD:
         return "INFEASIBLE or UNBOUNDED";
      case CPX_STAT_NUM_BEST:
         return "CPX_STAT_NUM_BEST: SOL AVAILABLE BUT NOT PROVEN OPTIMAL DUE TO NUM TROUBLE";
      default:
         return "UNKNOWN";
      }  /*lint !e788*/

      return "UNKNOWN";
   }
#endif

#ifndef NDEBUG
   bool checkConsistentBounds()
   {
      for( int i = 0; i < nCols(); ++i )
      {
         if( lower(i) > upper(i) )
         {
            SCIPerrorMessage("inconsistent bounds on column %d: lower=%.17g, upper=%.17g\n",
               i, lower(i), upper(i));
            return false;
         }
      }

      return true;
   }

   bool checkConsistentSides()
   {
      for( int i = 0; i < nRows(); ++i )
      {
         if( lhs(i) > rhs(i) )
         {
            SCIPerrorMessage("inconsistent sides on row %d: lhs=%.17g, rhs=%.17g\n",
               i, lhs(i), rhs(i));
            return false;
         }
      }

      return true;
   }
#endif

   void trySolve(bool printwarning = true)
   {
      Real timespent;
      Real timelimit;
      try
      {
	 m_stat = SPxSolver::solve();
      }
      catch(SPxException x)
      {
	 std::string s = x.what();
         if( printwarning )
         {
            SCIPmessagePrintWarning(m_messagehdlr, "SoPlex threw an exception: %s\n", s.c_str());
         }
	 m_stat = SPxSolver::status();

	 /* since it is not clear if the status in SoPlex are set correctly
	  * we want to make sure that if an error is thrown the status is
          * not OPTIMAL anymore.
          */
	 assert( m_stat != SPxSolver::OPTIMAL );
      }

      /* save iteration count */
      m_itused += SPxSolver::iterations();
      assert(m_itlim < 0 || m_itused <= m_itlim);

      /* update time limit */
      timespent = SPxSolver::time();
      if( timespent > 0 )
      {
         /* get current time limit */
         timelimit = SPxSolver::terminationTime();
         if( timelimit > timespent )
            timelimit -= timespent;
         else
            timelimit = 0;
         /* set new time limit */
         assert(timelimit >= 0);
         SPxSolver::setTerminationTime(timelimit);
      }
   }

   virtual Status doSolve(bool printwarning = true)
   {
      int verbosity;

      /* store and set verbosity */
      verbosity = Param::verbose();
      Param::setVerbose(getLpInfo() ? SOPLEX_VERBLEVEL : 0);

      assert(checkConsistentBounds());
      assert(checkConsistentSides());

#ifdef WITH_LPSCHECK
      /* dump LP with current basis and settings saved in SoPlex */
      if( getDoubleCheck() )
         writeState("spxcheck", NULL, NULL);
#endif

      /* in auto pricing, do the first iterations with devex, then switch to steepest edge */
      setTerminationIter(m_autopricing && (m_itlim < 0 || m_itlim - m_itused > AUTOPRICING_ITERSWITCH) ? AUTOPRICING_ITERSWITCH : m_itlim - m_itused);

      trySolve(printwarning);

      if( m_autopricing && m_stat == SPxSolver::ABORT_ITER && (m_itlim < 0 || m_itlim - m_itused > 0) )
      {
         setTerminationIter(m_itlim - m_itused);
         setPricer(&m_price_steep_ex);

         trySolve(printwarning);

         setPricer(&m_price_devex);
      }

      /* for safety reset iteration limit */
      setTerminationIter(m_itlim);

      if( m_stat == OPTIMAL )
      {
         Real objval = value();

         if( (objval > m_objUpLimit) || (objval < m_objLoLimit) )
            m_stat = ABORT_VALUE;
      }

#ifdef WITH_LPSCHECK
      /* if SoPlex gave a definite answer, we double check if it is consistent with CPLEX's answer */
      if( getDoubleCheck() && (m_stat == SPxSolver::OPTIMAL || m_stat == SPxSolver::UNBOUNDED || m_stat == SPxSolver::INFEASIBLE || m_stat == SPxSolver::ABORT_VALUE) )
      {
         SCIP_Real cpxobj;
         int cpxstat;

         /* read LP with basis */
         CPX_CALL( CPXreadcopyprob(m_cpxenv, m_cpxlp, "spxcheck.mps", NULL) );
         CPX_CALL( CPXreadcopybase(m_cpxenv, m_cpxlp, "spxcheck.bas") );

         /* set tolerances */
         CPX_CALL( CPXsetdblparam(m_cpxenv, CPX_PARAM_EPOPT, MAX(opttol(), 1e-9)) );
         CPX_CALL( CPXsetdblparam(m_cpxenv, CPX_PARAM_EPRHS, MAX(feastol(), 1e-9)) );

         /* solve LP */
         CPX_CALL( CPXlpopt(m_cpxenv, m_cpxlp) );

         /* get solution status and objective value */
         CPX_CALL( CPXsolution(m_cpxenv, m_cpxlp, &cpxstat, &cpxobj, NULL, NULL, NULL, NULL) );
         if( getSense() == SPxLP::MAXIMIZE )
            cpxobj *= -1.0;

         /* check for inconsistent statuses */
         if( cpxstat == CPX_STAT_OPTIMAL_INFEAS )
         {
            SCIPerrorMessage("In %s: SoPlex status=%d (%s) while CPLEX status=%d (%s)\n",
               m_probname, m_stat, spxStatusString(m_stat), cpxstat, cpxStatusString(cpxstat));
            if( EXIT_AT_CPXERROR )
               exit(1);
         }
         else if( (m_stat == SPxSolver::OPTIMAL && cpxstat != CPX_STAT_OPTIMAL)
            || (m_stat == SPxSolver::UNBOUNDED && cpxstat != CPX_STAT_UNBOUNDED)
            || (m_stat == SPxSolver::INFEASIBLE && cpxstat != CPX_STAT_INFEASIBLE) )
         {
            SCIPerrorMessage("In %s: SoPlex status=%d (%s) while CPLEX status=%d (%s) (checknum=%d)\n",
               m_probname, m_stat, spxStatusString(m_stat), cpxstat, cpxStatusString(cpxstat), m_checknum);
            if( EXIT_AT_WRONG_RESULT )
               exit(1);
         }
         else if( m_stat == SPxSolver::ABORT_VALUE )
         {
            switch( cpxstat )
            {
            case CPX_STAT_OPTIMAL:
               if( (getSense() == SPxSolver::MINIMIZE && LTrel(cpxobj, getObjUpLimit(), 2*opttol()))
                  || (getSense() == SPxSolver::MAXIMIZE && GTrel(cpxobj, getObjLoLimit(), 2*opttol())) )
               {
                  SCIPerrorMessage("In %s: SoPlex returned status=%d (%s) while CPLEX claims obj=%.10f %s %.10f=obj.limit (%s) (checknum=%d)\n",
                     m_probname, m_stat, spxStatusString(m_stat), cpxobj, getSense() == SPxSolver::MINIMIZE ? "<" : ">",
                     getSense() == SPxSolver::MINIMIZE ? getObjUpLimit() : getObjLoLimit(), cpxStatusString(cpxstat), m_checknum);
                  if( EXIT_AT_WRONG_RESULT )
                     exit(1);
               }
               else if( (getSense() == SPxSolver::MINIMIZE && cpxobj < getObjUpLimit())
                  || (getSense() == SPxSolver::MAXIMIZE && cpxobj > getObjLoLimit()) )
               {
                  SCIPerrorMessage("In %s: SoPlex returned status=%d (%s) while CPLEX claims obj=%.10f %s %.10f=obj.limit (%s) (checknum=%d)\n",
                     m_probname, m_stat, spxStatusString(m_stat), cpxobj, getSense() == SPxSolver::MINIMIZE ? "<" : ">",
                     getSense() == SPxSolver::MINIMIZE ? getObjUpLimit() : getObjLoLimit(), cpxStatusString(cpxstat), m_checknum);
               }
               break;
            case CPX_STAT_OPTIMAL_INFEAS:
            case CPX_STAT_NUM_BEST:
               if( (getSense() == SPxSolver::MINIMIZE && cpxobj < getObjUpLimit())
                  || (getSense() == SPxSolver::MAXIMIZE && cpxobj > getObjLoLimit()) )
               {
                  SCIPerrorMessage("In %s: SoPlex returned status=%d (%s) while CPLEX claims obj=%.10f %s %.10f=obj.limit (%s) (checknum=%d)\n",
                     m_probname, m_stat, spxStatusString(m_stat), cpxobj, getSense() == SPxSolver::MINIMIZE ? "<" : ">",
                     getSense() == SPxSolver::MINIMIZE ? getObjUpLimit() : getObjLoLimit(), cpxStatusString(cpxstat), m_checknum);
               }
               break;
            case CPX_STAT_INFEASIBLE:
               break;
            case CPX_STAT_UNBOUNDED:
               SCIPerrorMessage("In %s: SoPlex status=%d (%s) while CPLEX status=%d (%s) (checknum=%d)\n",
                  m_probname, m_stat, spxStatusString(m_stat), cpxstat, cpxStatusString(cpxstat), m_checknum);
               if( EXIT_AT_WRONG_RESULT )
                  exit(1);
               break;
            case CPX_STAT_INForUNBD:
            default:
               SCIPerrorMessage("In %s: SoPlex status=%d (%s) while CPLEX status=%d (%s) (checknum=%d)\n",
                  m_probname, m_stat, spxStatusString(m_stat), cpxstat, cpxStatusString(cpxstat), m_checknum);
               break;
            }  /*lint !e788*/
         }
         /* check for same objective values */
         else if( m_stat == SPxSolver::OPTIMAL )
         {
            if( (getSense() == SPxSolver::MINIMIZE && LTrel(value(), cpxobj, 2*opttol()))
               || (getSense() == SPxSolver::MAXIMIZE && GTrel(value(), cpxobj, 2*opttol())) )
            {
               /* SCIPerrorMessage("In %s: LP optimal; SoPlex value=%.10f %s CPLEX value=%.10f too good (checknum=%d)\n", value(),
                  m_probname, getSense() == SPxSolver::MINIMIZE ? "<" : ">", cpxobj, m_checknum); */
            }
            else if( (getSense() == SPxSolver::MINIMIZE && GTrel(value(), cpxobj, 2*opttol()))
               || (getSense() == SPxSolver::MAXIMIZE && LTrel(value(), cpxobj, 2*opttol())) )
            {
               SCIPerrorMessage("In %s: LP optimal; SoPlex value=%.10f %s CPLEX value=%.10f suboptimal (checknum=%d)\n", value(),
                  m_probname, getSense() == SPxSolver::MINIMIZE ? ">" : "<", cpxobj, m_checknum);
               if( EXIT_AT_WRONG_RESULT )
                  exit(1);
            }
         }
      }

   ENDCHECK:
#endif

      /* restore verbosity */
      Param::setVerbose(verbosity);

      return m_stat;
   }

   virtual Status solve()
   {
      assert(m_sense == sense());

      SPxEquiliSC* scaler = NULL;
      SPxSimplifier* simplifier = NULL;
      SPxLP origlp;
      SPxSimplifier::Result result = SPxSimplifier::OKAY;

      /* delete starting basis if solving from scratch */
      if ( getFromScratch() )
      {
	 try
	 {
	    SPxSolver::reLoad();
	 }
	 catch(SPxException x)
	 {
	    std::string s = x.what();
	    SCIPmessagePrintWarning(m_messagehdlr, "SoPlex threw an exception: %s\n", s.c_str());
	    m_stat = SPxSolver::status();
	    assert( m_stat != SPxSolver::OPTIMAL );
	    return m_stat;
	 }
      }
      assert(!getFromScratch() || getBasisStatus() == SPxBasis::NO_PROBLEM);

      /* use scaler and simplifier if no basis is loaded, i.e., if solving for the first time or from scratch */
      if( SPxSolver::getBasisStatus() == SPxBasis::NO_PROBLEM && getScaling() && nCols() > 0 && nRows() > 0 )
      {
         scaler = new SPxEquiliSC();
         assert(scaler != NULL);
      }

      if( SPxSolver::getBasisStatus() == SPxBasis::NO_PROBLEM && getPresolving() && nCols() > 0 && nRows() > 0 )
      {
         simplifier = new SPxMainSM();
         assert(simplifier != NULL);
      }

      /* store original lp */
      if( scaler != NULL || simplifier != NULL )
         origlp = SPxLP(*this);

   SOLVEAGAIN:
      /* perform scaling and presolving */
      if( scaler != NULL )
      {
         SCIPdebugMessage("scaling LP\n");
         scaler->scale(*this);
      }

      if( simplifier != NULL )
      {
         int verbosity;

         /* store and set verbosity */
         verbosity = Param::verbose();
         Param::setVerbose(getLpInfo() ? SOPLEX_VERBLEVEL : 0);
         SCIPdebugMessage("simplifying LP\n");
#if ((SOPLEX_VERSION == 160 && SOPLEX_SUBVERSION >= 5) || SOPLEX_VERSION > 160)
         result = simplifier->simplify(*this, epsilon(), feastol(), opttol());
#else
         result = simplifier->simplify(*this, epsilon(), delta());
#endif
         SCIPdebugMessage("simplifier ended with status %u (0: OKAY, 1: INFEASIBLE, 2: DUAL_INFEASIBLE, 3: UNBOUNDED, 4: VANISHED)\n", result);

         /* unsimplification is not designed for these cases, thus reload original/scaled lp */
         if( result == SPxSimplifier::INFEASIBLE || result == SPxSimplifier::DUAL_INFEASIBLE )
         {
            SCIPdebugMessage("simplifier detected primal or dual infeasibility - reloading and solving unsimplified LP\n");

            delete simplifier;
            simplifier = NULL;

            SPxSolver::loadLP(origlp);
            m_sense = sense();

            goto SOLVEAGAIN;
         }
         /* reset verbosity */
         Param::setVerbose(verbosity);
      }

      /* solve */
      m_itused = 0;
      if( result != SPxSimplifier::VANISHED )
      {
         /* we have to deactivate the objective limit, since we do not know the transformed value */
         Real objlolimit = getObjLoLimit();
         Real objuplimit = getObjUpLimit();

         if( simplifier != NULL || scaler != NULL )
         {
            setObjLoLimit(-soplex::infinity);
            setObjUpLimit(soplex::infinity);
         }

#ifndef NDEBUG
         doSolve();
#else
         doSolve(false);
#endif

         if( simplifier != NULL || scaler != NULL )
         {
            setObjLoLimit(objlolimit);
            setObjUpLimit(objuplimit);
         }
      }

      /* unsimplification only stable for optimal basis */
      if( m_stat != SPxSolver::OPTIMAL && simplifier != NULL )
      {
         SCIPdebugMessage("presolved LP not optimal - reloading and solving original LP\n");

         delete simplifier;
         simplifier = NULL;

         SPxSolver::loadLP(origlp);
         m_sense = sense();

         goto SOLVEAGAIN;
      }

      /* if scaling or presolving was applied, restore original lp */
      if( scaler != NULL || simplifier != NULL )
      {
         SPxSolver::VarStatus* cstat = NULL;
         SPxSolver::VarStatus* rstat = NULL;

         /* get basis if at least regular */
         if( (simplifier == NULL || result != SPxSimplifier::VANISHED) && SPxSolver::getBasisStatus() >= SPxBasis::REGULAR )
         {
            SCIPdebugMessage("get basis of presolved LP\n");
            rstat = new SPxSolver::VarStatus[nRows()];
            cstat = new SPxSolver::VarStatus[nCols()];
            SPxSolver::getBasis(rstat, cstat);
         }

         /* unsimplify */
         if( simplifier != NULL && SPxSolver::getBasisStatus() >= SPxBasis::REGULAR )
         {
            assert((result == SPxSimplifier::VANISHED) == (cstat == NULL));
            assert((result == SPxSimplifier::VANISHED) == (rstat == NULL));

            /* dimension of presolved lp */
            int ncols = result == SPxSimplifier::VANISHED ? 0 : nCols();
            int nrows = result == SPxSimplifier::VANISHED ? 0 : nRows();

            /* get solution of presolved lp */
            DVector primals(ncols);
            DVector duals(nrows);
            DVector slacks(nrows);
            DVector redcosts(ncols);
            if( result != SPxSimplifier::VANISHED )
            {
               SPxSolver::getPrimal(primals);
               SPxSolver::getDual(duals);
               SPxSolver::getSlacks(slacks);
               SPxSolver::getRedCost(redcosts);
            }

            /* perform unsimplification */
            SCIPdebugMessage("unsimplify\n");
            try
            {
               simplifier->unsimplify(primals, duals, slacks, redcosts, rstat, cstat);
            }
            catch(SPxException x)
            {
               std::string s = x.what();
               SCIPmessagePrintWarning(m_messagehdlr, "SoPlex unsimplification unsuccessful; solving again without LP presolving (SoPlex says %s)\n",
                  s.c_str());
            }

            if( cstat != NULL )
            {
               delete[] cstat;
               cstat = NULL;
            }
            if( rstat != NULL )
            {
               delete[] rstat;
               rstat = NULL;
            }

            if( simplifier->isUnsimplified() )
            {
               /* get basis for original lp */
               rstat = new SPxSolver::VarStatus[origlp.nRows()];
               cstat = new SPxSolver::VarStatus[origlp.nCols()];
               simplifier->getBasis(rstat, cstat);
            }
         }

         /* reload original lp */
         SCIPdebugMessage("reload original LP\n");
         SPxSolver::loadLP(origlp);
         m_sense = sense();

         /* set basis from preprocessed lp and reoptimize */
         if( rstat != NULL && cstat != NULL )
         {
            SCIPdebugMessage("load unsimplified basis into original LP\n");
            SPxSolver::setBasis(rstat, cstat);
         }

         SCIPdebugMessage("solve original LP\n");
#ifndef NDEBUG
         doSolve();
#else
         doSolve(false);
#endif

         /* free allocated memory */
         if( cstat != NULL )
            delete[] cstat;
         if( rstat != NULL )
            delete[] rstat;
         if( scaler != NULL )
            delete scaler;
         if( simplifier != NULL )
            delete simplifier;
      }

      if( m_stat == OPTIMAL )
      {
         Real objval = value();
         
         if( (objval > m_objUpLimit) || (objval < m_objLoLimit) )
            m_stat = ABORT_VALUE;
      }

      return m_stat;
   }

   /** save the current basis */
   void savePreStrongbranchingBasis()
   {
      assert(m_rowstat == NULL);
      assert(m_colstat == NULL);

      m_rowstat = new SPxSolver::VarStatus[nRows()];
      m_colstat = new SPxSolver::VarStatus[nCols()];

      try
      {
         m_stat = getBasis(m_rowstat, m_colstat);
      }
      catch(SPxException x)
      {
#ifndef NDEBUG
         std::string s = x.what();
         SCIPmessagePrintWarning(m_messagehdlr, "SoPlex threw an exception: %s\n", s.c_str());

         /* since it is not clear if the status in SoPlex are set correctly
          * we want to make sure that if an error is thrown the status is
          * not OPTIMAL anymore.
          */
         assert(m_stat != SPxSolver::OPTIMAL);
#endif
      }
   }

   /** restore basis */
   void restorePreStrongbranchingBasis()
   {
      assert(m_rowstat != NULL);
      assert(m_colstat != NULL);

      try
      {
         setBasis(m_rowstat, m_colstat);
      }
      catch(SPxException x)
      {
#ifndef NDEBUG
         std::string s = x.what();
         SCIPmessagePrintWarning(m_messagehdlr, "SoPlex threw an exception: %s\n", s.c_str());
#endif
         m_stat = SPxSolver::status();

         /* since it is not clear if the status in SoPlex are set correctly
          * we want to make sure that if an error is thrown the status is
          * not OPTIMAL anymore.
          */
         assert(m_stat != SPxSolver::OPTIMAL);
      }
   }

   /** if basis is in store, delete it without restoring it */
   void freePreStrongbranchingBasis()
   {
      if ( m_rowstat != NULL )
      {
         delete [] m_rowstat;
         m_rowstat = NULL;
      }
      if ( m_colstat != NULL ) 
      {
         delete [] m_colstat;
         m_colstat = NULL;
      }
   }

   /** is pre-strong-branching basis freed? */
   bool preStrongbranchingBasisFreed()
   {
      return ((m_rowstat == NULL ) && (m_colstat == NULL));
   }

   Status getStatus() const
   {
      return m_stat;
   }

   Status updateStatus()
   {
      m_stat = SPxSolver::status();
      return m_stat;
   }

   bool isInitialized() const
   {
      return SPxSolver::isInitialized();
   }

   int iterations() const
   {
      return m_itused;
   }

   virtual void clear()
   {
      SPxSolver::clear();
      freePreStrongbranchingBasis();
      m_stat = NO_PROBLEM;
      m_sense = sense();
   }

   bool readLP(const char* fname)
   {
      clear();

      if ( m_rownames != 0 )
         delete m_rownames;
      if ( m_colnames != 0 )
         delete m_colnames;
      m_rownames = new NameSet;
      m_colnames = new NameSet;

      if( SPxSolver::readFile(fname, m_rownames, m_colnames) )
      {
         m_stat = NO_PROBLEM;
         m_sense = sense();
         return true;
      }

      return false;
   }

   /** copy column names into namestorage with access via colnames */
   void getColNames(
      int                firstcol,           /**< first column to get name from LP */
      int                lastcol,            /**< last column to get name from LP */
      char**             colnames,           /**< pointers to column names (of size at least lastcol-firstcol+1) */
      char*              namestorage,        /**< storage for col names */
      int                namestoragesize,    /**< size of namestorage (if 0, storageleft returns the storage needed) */
      int*               storageleft         /**< amount of storage left (if < 0 the namestorage was not big enough) */
      )
   {
      assert( m_colnames != NULL );

      // compute size
      if ( namestoragesize == 0 )
      {
	 // the following may overestimate the space requirements
	 *storageleft = -m_colnames->memSize();
      }
      else
      {
	 NameSet* names = m_colnames;
	 assert( names != 0 );
	 int sizeleft = namestoragesize;
	 char* s = namestorage;
	 for (int j = firstcol; j <= lastcol; ++j)
	 {
	    const char* t = (*names)[j];
	    colnames[j-firstcol] = s;
	    while( *t != '\0' && sizeleft >= 0 )
	    {
	       *(s++) = *(t++);
	       --sizeleft;
	    }
	    *(s++) = '\0';
	 }
	 if ( sizeleft == 0 )
	 {
	    *storageleft = namestoragesize - m_colnames->memSize();
	    assert( *storageleft <= 0 );
	 }
	 else
	    *storageleft = sizeleft;
      }
   }

   /** copy row names into namestorage with access via row */
   void getRowNames(
      int                firstrow,           /**< first row to get name from LP */
      int                lastrow,            /**< last row to get name from LP */
      char**             rownames,           /**< pointers to row names (of size at least lastrow-firstrow+1) */
      char*              namestorage,        /**< storage for row names */
      int                namestoragesize,    /**< size of namestorage (if 0, -storageleft returns the storage needed) */
      int*               storageleft         /**< amount of storage left (if < 0 the namestorage was not big enough) */
      )
   {
      assert( m_rownames != NULL );

      // compute size
      if ( namestoragesize == 0 )
      {
	 // the following may overestimate the space requirements
	 *storageleft = -m_rownames->memSize();
      }
      else
      {
	 NameSet* names = m_rownames;
	 assert( names != 0 );
	 int sizeleft = namestoragesize;
	 char* s = namestorage;
	 for (int i = firstrow; i <= lastrow; ++i)
	 {
	    const char* t = (*names)[i];
	    rownames[i-firstrow] = s;
	    while( *t != '\0' && sizeleft >= 0 )
	    {
	       *(s++) = *(t++);
	       --sizeleft;
	    }
	    *(s++) = '\0';
	 }
	 if ( sizeleft == 0 )
	 {
	    *storageleft = m_rownames->memSize() - namestoragesize;
	    assert( *storageleft <= 0 );
	 }
	 else
	    *storageleft = sizeleft;
      }
   }
}; /*lint !e1748*/




/********************************************************************/
/*-----------------------------  C  --------------------------------*/
/********************************************************************/

#include "scip/lpi.h"
#include "scip/bitencode.h"

typedef SCIP_DUALPACKET COLPACKET;           /* each column needs two bits of information (basic/on_lower/on_upper) */
#define COLS_PER_PACKET SCIP_DUALPACKETSIZE
typedef SCIP_DUALPACKET ROWPACKET;           /* each row needs two bit of information (basic/on_lower/on_upper) */
#define ROWS_PER_PACKET SCIP_DUALPACKETSIZE



/** LP interface */
struct SCIP_LPi
{
   SPxSCIP*              spx;                /**< our SPxSolver implementation */
   int*                  cstat;              /**< array for storing column basis status */
   int*                  rstat;              /**< array for storing row basis status */
   int                   cstatsize;          /**< size of cstat array */
   int                   rstatsize;          /**< size of rstat array */
   SCIP_PRICING          pricing;            /**< current pricing strategy */
   SCIP_Bool             solved;             /**< was the current LP solved? */
   SLUFactor*            factorization;      /**< factorization possibly needed for basis inverse */
   SCIP_Real             rowrepswitch;       /**< use row representation if number of rows divided by number of columns exceeds this value */
   SCIP_MESSAGEHDLR*     messagehdlr;        /**< messagehdlr handler to printing messages, or NULL */
};

/** LPi state stores basis information */
struct SCIP_LPiState
{
   int                   ncols;              /**< number of LP columns */
   int                   nrows;              /**< number of LP rows */
   COLPACKET*            packcstat;          /**< column basis status in compressed form */
   ROWPACKET*            packrstat;          /**< row basis status in compressed form */
};




/*
 * dynamic memory arrays
 */

/** resizes cstat array to have at least num entries */
static
SCIP_RETCODE ensureCstatMem(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   int                   num                 /**< minimal number of entries in array */
   )
{
   assert(lpi != NULL);

   if( num > lpi->cstatsize )
   {
      int newsize;

      newsize = MAX(2*lpi->cstatsize, num);
      SCIP_ALLOC( BMSreallocMemoryArray(&lpi->cstat, newsize) );
      lpi->cstatsize = newsize;
   }
   assert(num <= lpi->cstatsize);

   return SCIP_OKAY;
}

/** resizes rstat array to have at least num entries */
static
SCIP_RETCODE ensureRstatMem(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   int                   num                 /**< minimal number of entries in array */
   )
{
   assert(lpi != NULL);

   if( num > lpi->rstatsize )
   {
      int newsize;

      newsize = MAX(2*lpi->rstatsize, num);
      SCIP_ALLOC( BMSreallocMemoryArray(&lpi->rstat, newsize) );
      lpi->rstatsize = newsize;
   }
   assert(num <= lpi->rstatsize);

   return SCIP_OKAY;
}




/*
 * LPi state methods
 */

/** returns the number of packets needed to store column packet information */
static 
int colpacketNum(
   int                   ncols               /**< number of columns to store */
   )
{
   return (ncols+(int)COLS_PER_PACKET-1)/(int)COLS_PER_PACKET;
}

/** returns the number of packets needed to store row packet information */
static 
int rowpacketNum(
   int                   nrows               /**< number of rows to store */
   )
{
   return (nrows+(int)ROWS_PER_PACKET-1)/(int)ROWS_PER_PACKET;
}

/** store row and column basis status in a packed LPi state object */
static
void lpistatePack(
   SCIP_LPISTATE*       lpistate,            /**< pointer to LPi state data */
   const int*           cstat,               /**< basis status of columns in unpacked format */
   const int*           rstat                /**< basis status of rows in unpacked format */
   )
{
   assert(lpistate != NULL);
   assert(lpistate->packcstat != NULL);
   assert(lpistate->packrstat != NULL);

   SCIPencodeDualBit(cstat, lpistate->packcstat, lpistate->ncols);
   SCIPencodeDualBit(rstat, lpistate->packrstat, lpistate->nrows);
}

/** unpacks row and column basis status from a packed LPi state object */
static
void lpistateUnpack(
   const SCIP_LPISTATE* lpistate,            /**< pointer to LPi state data */
   int*                 cstat,               /**< buffer for storing basis status of columns in unpacked format */
   int*                 rstat                /**< buffer for storing basis status of rows in unpacked format */
   )
{
   assert(lpistate != NULL);
   assert(lpistate->packcstat != NULL);
   assert(lpistate->packrstat != NULL);

   SCIPdecodeDualBit(lpistate->packcstat, cstat, lpistate->ncols);
   SCIPdecodeDualBit(lpistate->packrstat, rstat, lpistate->nrows);
}

/** creates LPi state information object */
static
SCIP_RETCODE lpistateCreate(
   SCIP_LPISTATE**       lpistate,           /**< pointer to LPi state */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   int                   ncols,              /**< number of columns to store */
   int                   nrows               /**< number of rows to store */
   )
{
   assert(lpistate != NULL);
   assert(blkmem != NULL);
   assert(ncols >= 0);
   assert(nrows >= 0);

   SCIP_ALLOC( BMSallocBlockMemory(blkmem, lpistate) );
   SCIP_ALLOC( BMSallocBlockMemoryArray(blkmem, &(*lpistate)->packcstat, colpacketNum(ncols)) );
   SCIP_ALLOC( BMSallocBlockMemoryArray(blkmem, &(*lpistate)->packrstat, rowpacketNum(nrows)) );

   return SCIP_OKAY;
}

/** frees LPi state information */
static
void lpistateFree(
   SCIP_LPISTATE**       lpistate,           /**< pointer to LPi state information (like basis information) */
   BMS_BLKMEM*           blkmem              /**< block memory */
   )
{
   assert(blkmem != NULL);
   assert(lpistate != NULL);
   assert(*lpistate != NULL);

   BMSfreeBlockMemoryArray(blkmem, &(*lpistate)->packcstat, colpacketNum((*lpistate)->ncols));
   BMSfreeBlockMemoryArray(blkmem, &(*lpistate)->packrstat, rowpacketNum((*lpistate)->nrows));
   BMSfreeBlockMemory(blkmem, lpistate);
}




/*
 * local methods
 */

/** converts SCIP's objective sense into SoPlex's objective sense */
static
SPxLP::SPxSense spxObjsen(
   SCIP_OBJSEN           objsen              /**< SCIP's objective sense value */
   )
{
   switch( objsen )
   {
   case SCIP_OBJSEN_MAXIMIZE:
      return SPxLP::MAXIMIZE;
   case SCIP_OBJSEN_MINIMIZE:
      return SPxLP::MINIMIZE;
   default:
      SCIPerrorMessage("invalid objective sense\n");
      SCIPABORT();
      return SPxLP::MINIMIZE;
   }
}

/** marks the current LP to be unsolved */
static
void invalidateSolution(SCIP_LPI* lpi)
{
   assert(lpi != NULL);
   lpi->solved = FALSE;
   if ( lpi->factorization != 0 )
   {
      delete lpi->factorization;
      lpi->factorization = 0;
   }
}



/*
 * LP Interface Methods
 */


/*
 * Miscellaneous Methods
 */

static char spxname[100];
static char spxdesc[200];

/**@name Miscellaneous Methods */
/**@{ */

/** gets name and version of LP solver */
const char* SCIPlpiGetSolverName(
   void
   )
{
   SCIPdebugMessage("calling SCIPlpiGetSolverName()\n");

#if (SOPLEX_SUBVERSION > 0)
   sprintf(spxname, "SoPlex %d.%d.%d.%d", SOPLEX_VERSION/100, (SOPLEX_VERSION % 100)/10, SOPLEX_VERSION % 10, SOPLEX_SUBVERSION);
#else
   sprintf(spxname, "SoPlex %d.%d.%d", SOPLEX_VERSION/100, (SOPLEX_VERSION % 100)/10, SOPLEX_VERSION % 10);
#endif
   return spxname;
}

/** gets description of LP solver (developer, webpage, ...) */
const char* SCIPlpiGetSolverDesc(
   void
   )
{
   sprintf(spxdesc, "%s", "Linear Programming Solver developed at Zuse Institute Berlin (soplex.zib.de)");
#if (SOPLEX_VERSION >= 160)
   sprintf(spxdesc, "%s [GitHash: %s]", spxdesc, getGitHash());
#endif
#ifdef WITH_LPSCHECK
   sprintf(spxdesc, "%s %s", spxdesc, "- including CPLEX double check");
#endif
   return spxdesc;
}

/** gets pointer for LP solver - use only with great care */
void* SCIPlpiGetSolverPointer(
   SCIP_LPI*             lpi                 /**< pointer to an LP interface structure */
   )
{
   return (void*) lpi->spx;
}
/**@} */




/*
 * LPI Creation and Destruction Methods
 */

/**@name LPI Creation and Destruction Methods */
/**@{ */

/** creates an LP problem object */
SCIP_RETCODE SCIPlpiCreate(
   SCIP_LPI**            lpi,                /**< pointer to an LP interface structure */
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler to use for printing messages, or NULL */
   const char*           name,               /**< problem name */
   SCIP_OBJSEN           objsen              /**< objective sense */
   )
{
   assert(lpi != NULL);

   /* create SoPlex object */
   SCIP_ALLOC( BMSallocMemory(lpi) );

   /* we use this construction to allocate the memory for the SoPlex class also via the blockmemshell */
   (*lpi)->spx = static_cast<SPxSCIP*>(BMSallocMemoryCPP(sizeof(SPxSCIP)));
   SOPLEX_TRY( messagehdlr, (*lpi)->spx = new ((*lpi)->spx) SPxSCIP(messagehdlr, name) );
   (*lpi)->cstat = NULL;
   (*lpi)->rstat = NULL;
   (*lpi)->cstatsize = 0;
   (*lpi)->rstatsize = 0;
   (*lpi)->pricing = SCIP_PRICING_LPIDEFAULT;
   (*lpi)->factorization = 0;
   (*lpi)->rowrepswitch = SCIPlpiInfinity(*lpi);
   (*lpi)->messagehdlr = messagehdlr;

   invalidateSolution(*lpi);

   /* set objective sense */
   SCIP_CALL( SCIPlpiChgObjsen(*lpi, objsen) );

   /* set default pricing */
   SCIP_CALL( SCIPlpiSetIntpar(*lpi, SCIP_LPPAR_PRICING, (int)(*lpi)->pricing) );

   return SCIP_OKAY;
}

/** deletes an LP problem object */
SCIP_RETCODE SCIPlpiFree(
   SCIP_LPI**            lpi                 /**< pointer to an LP interface structure */
   )
{
   assert(lpi != NULL);
   assert(*lpi != NULL);
   assert((*lpi)->spx != NULL);

   /* free LP using destructor and free memory via blockmemshell */
   (*lpi)->spx->~SPxSCIP();
   BMSfreeMemory(&((*lpi)->spx));

   /* free memory */
   BMSfreeMemoryArrayNull(&(*lpi)->cstat);
   BMSfreeMemoryArrayNull(&(*lpi)->rstat);
   BMSfreeMemory(lpi);

   return SCIP_OKAY;
}

/**@} */




/*
 * Modification Methods
 */

/**@name Modification Methods */
/**@{ */

/** copies LP data with column matrix into LP solver */
SCIP_RETCODE SCIPlpiLoadColLP(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   SCIP_OBJSEN           objsen,             /**< objective sense */
   int                   ncols,              /**< number of columns */
   const SCIP_Real*      obj,                /**< objective function values of columns */
   const SCIP_Real*      lb,                 /**< lower bounds of columns */
   const SCIP_Real*      ub,                 /**< upper bounds of columns */
   char**                colnames,           /**< column names, or NULL */
   int                   nrows,              /**< number of rows */
   const SCIP_Real*      lhs,                /**< left hand sides of rows */
   const SCIP_Real*      rhs,                /**< right hand sides of rows */
   char**                /*rownames*/,       /**< row names, or NULL */
   int                   nnonz,              /**< number of nonzero elements in the constraint matrix */
   const int*            beg,                /**< start index of each column in ind- and val-array */
   const int*            ind,                /**< row indices of constraint matrix entries */
   const SCIP_Real*      val                 /**< values of constraint matrix entries */
   )
{
   SCIPdebugMessage("calling SCIPlpiLoadColLP()\n");

   assert(lpi != NULL);
   assert(lpi->spx != NULL);
   assert(lhs != NULL);
   assert(rhs != NULL);

   invalidateSolution(lpi);
   assert( lpi->spx->preStrongbranchingBasisFreed() );

   try
   {
      SPxSCIP* spx = lpi->spx;
      LPRowSet rows(nrows);
      DSVector emptyVector(0);
      int i;

      spx->clear();

      /* set objective sense */
      spx->setSense(spxObjsen(objsen));

      /* create empty rows with given sides */
      for( i = 0; i < nrows; ++i )
         rows.add(lhs[i], emptyVector, rhs[i]);
      spx->addRows(rows);
   
      /* create column vectors with coefficients and bounds */
      SCIP_CALL( SCIPlpiAddCols(lpi, ncols, obj, lb, ub, colnames, nnonz, beg, ind, val) );
   }
   catch(SPxException x)
   {
#ifndef NDEBUG
      std::string s = x.what();
      SCIPmessagePrintWarning(lpi->messagehdlr, "SoPlex threw an exception: %s\n", s.c_str());
#endif
      return SCIP_LPERROR;
   }

   return SCIP_OKAY;
}

/** adds columns to the LP */
SCIP_RETCODE SCIPlpiAddCols(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   int                   ncols,              /**< number of columns to be added */
   const SCIP_Real*      obj,                /**< objective function values of new columns */
   const SCIP_Real*      lb,                 /**< lower bounds of new columns */
   const SCIP_Real*      ub,                 /**< upper bounds of new columns */
   char**                /*colnames*/,       /**< column names, or NULL */
   int                   nnonz,              /**< number of nonzero elements to be added to the constraint matrix */
   const int*            beg,                /**< start index of each column in ind- and val-array, or NULL if nnonz == 0 */
   const int*            ind,                /**< row indices of constraint matrix entries, or NULL if nnonz == 0 */
   const SCIP_Real*      val                 /**< values of constraint matrix entries, or NULL if nnonz == 0 */
   )
{
   SCIPdebugMessage("calling SCIPlpiAddCols()\n");

   assert(lpi != NULL);
   assert(lpi->spx != NULL);
   assert(obj != NULL);
   assert(lb != NULL);
   assert(ub != NULL);
   assert(nnonz == 0 || beg != NULL);
   assert(nnonz == 0 || ind != NULL);
   assert(nnonz == 0 || val != NULL);

   invalidateSolution(lpi);

   assert( lpi->spx->preStrongbranchingBasisFreed() );

   SPxSCIP* spx = lpi->spx;
   try
   {
      LPColSet cols(ncols);
      DSVector colVector(ncols);
      int start;
      int last;
      int i;

      /* create column vectors with coefficients and bounds */
      for( i = 0; i < ncols; ++i )
      {
         colVector.clear();
         if( nnonz > 0 )
         {
            start = beg[i];
            last = (i == ncols-1 ? nnonz : beg[i+1]);
            colVector.add( last-start, &ind[start], &val[start] );
         }
         cols.add(obj[i], lb[i], colVector, ub[i]);
      }
      spx->addCols(cols);
   }
   catch(SPxException x)
   {
#ifndef NDEBUG
      std::string s = x.what();      
      SCIPmessagePrintWarning(lpi->messagehdlr, "SoPlex threw an exception: %s\n", s.c_str());
#endif
      return SCIP_LPERROR;
   }
 
   return SCIP_OKAY;
}

/** deletes all columns in the given range from LP */
SCIP_RETCODE SCIPlpiDelCols(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   int                   firstcol,           /**< first column to be deleted */
   int                   lastcol             /**< last column to be deleted */
   )
{
   SCIPdebugMessage("calling SCIPlpiDelCols()\n");

   assert(lpi != NULL);
   assert(lpi->spx != NULL);
   assert(0 <= firstcol && firstcol <= lastcol && lastcol < lpi->spx->nCols());

   invalidateSolution(lpi);

   assert( lpi->spx->preStrongbranchingBasisFreed() );

   SOPLEX_TRY( lpi->messagehdlr, lpi->spx->removeColRange(firstcol, lastcol) );

   return SCIP_OKAY;
}

/** deletes columns from SCIP_LP; the new position of a column must not be greater that its old position */
SCIP_RETCODE SCIPlpiDelColset(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   int*                  dstat               /**< deletion status of columns
                                              *   input:  1 if column should be deleted, 0 if not
                                              *   output: new position of column, -1 if column was deleted */
   )
{
   int ncols;
   int i;

   SCIPdebugMessage("calling SCIPlpiDelColset()\n");

   assert(lpi != NULL);
   assert(lpi->spx != NULL);

   invalidateSolution(lpi);

   assert( lpi->spx->preStrongbranchingBasisFreed() );

   ncols = lpi->spx->nCols();

   /* SoPlex removeCols() method deletes the columns with dstat[i] < 0, so we have to negate the values */
   for( i = 0; i < ncols; ++i )
      dstat[i] *= -1;

   SOPLEX_TRY( lpi->messagehdlr, lpi->spx->removeCols(dstat) );

   return SCIP_OKAY;
}

/** adds rows to the LP */
SCIP_RETCODE SCIPlpiAddRows(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   int                   nrows,              /**< number of rows to be added */
   const SCIP_Real*      lhs,                /**< left hand sides of new rows */
   const SCIP_Real*      rhs,                /**< right hand sides of new rows */
   char**                /*rownames*/,       /**< row names, or NULL */
   int                   nnonz,              /**< number of nonzero elements to be added to the constraint matrix */
   const int*            beg,                /**< start index of each row in ind- and val-array, or NULL if nnonz == 0 */
   const int*            ind,                /**< column indices of constraint matrix entries, or NULL if nnonz == 0 */
   const SCIP_Real*      val                 /**< values of constraint matrix entries, or NULL if nnonz == 0 */
   )
{
   SCIPdebugMessage("calling SCIPlpiAddRows()\n");

   assert(lpi != NULL);
   assert(lpi->spx != NULL);
   assert(lhs != NULL);
   assert(rhs != NULL);
   assert(nnonz == 0 || beg != NULL);
   assert(nnonz == 0 || ind != NULL);
   assert(nnonz == 0 || val != NULL);

   invalidateSolution(lpi);

   assert( lpi->spx->preStrongbranchingBasisFreed() );

   try
   {
      SPxSCIP* spx = lpi->spx;
      LPRowSet rows(nrows);
      DSVector rowVector;
      int start;
      int last;
      int i;

      /* create row vectors with given sides */
      for( i = 0; i < nrows; ++i )
      {
         rowVector.clear();
         if( nnonz > 0 )
         {
            start = beg[i];
            last = (i == nrows-1 ? nnonz : beg[i+1]);
            rowVector.add( last-start, &ind[start], &val[start] );
         }
         rows.add(lhs[i], rowVector, rhs[i]);
      }
      spx->addRows(rows);
   }
   catch(SPxException x)
   {
#ifndef NDEBUG
      std::string s = x.what();      
      SCIPmessagePrintWarning(lpi->messagehdlr, "SoPlex threw an exception: %s\n", s.c_str());
#endif
      return SCIP_LPERROR;
   }

   return SCIP_OKAY;
}

/** deletes all rows in the given range from LP */
SCIP_RETCODE SCIPlpiDelRows(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   int                   firstrow,           /**< first row to be deleted */
   int                   lastrow             /**< last row to be deleted */
   )
{
   SCIPdebugMessage("calling SCIPlpiDelRows()\n");

   assert(lpi != NULL);
   assert(lpi->spx != NULL);
   assert(0 <= firstrow && firstrow <= lastrow && lastrow < lpi->spx->nRows());

   invalidateSolution(lpi);

   assert( lpi->spx->preStrongbranchingBasisFreed() );

   SOPLEX_TRY( lpi->messagehdlr, lpi->spx->removeRowRange(firstrow, lastrow) );

   return SCIP_OKAY;
}

/** deletes rows from SCIP_LP; the new position of a row must not be greater that its old position */
SCIP_RETCODE SCIPlpiDelRowset(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   int*                  dstat               /**< deletion status of rows
                                              *   input:  1 if row should be deleted, 0 if not
                                              *   output: new position of row, -1 if row was deleted */
   )
{
   int nrows;
   int i;

   SCIPdebugMessage("calling SCIPlpiDelRowset()\n");

   assert(lpi != NULL);
   assert(lpi->spx != NULL);

   invalidateSolution(lpi);

   assert( lpi->spx->preStrongbranchingBasisFreed() );

   nrows = lpi->spx->nRows();

   /* SoPlex removeRows() method deletes the rows with dstat[i] < 0, so we have to negate the values */
   for( i = 0; i < nrows; ++i )
      dstat[i] *= -1;

   SOPLEX_TRY( lpi->messagehdlr, lpi->spx->removeRows(dstat) );

   return SCIP_OKAY;
}

/** clears the whole LP */
SCIP_RETCODE SCIPlpiClear(
   SCIP_LPI*             lpi                 /**< LP interface structure */
   )
{
   SCIPdebugMessage("calling SCIPlpiClear()\n");

   assert(lpi != NULL);
   assert(lpi->spx != NULL);

   invalidateSolution(lpi);

   assert( lpi->spx->preStrongbranchingBasisFreed() );
   SOPLEX_TRY( lpi->messagehdlr, lpi->spx->clear() );

   return SCIP_OKAY;
}

/** changes lower and upper bounds of columns */
SCIP_RETCODE SCIPlpiChgBounds(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   int                   ncols,              /**< number of columns to change bounds for */
   const int*            ind,                /**< column indices */
   const SCIP_Real*      lb,                 /**< values for the new lower bounds */
   const SCIP_Real*      ub                  /**< values for the new upper bounds */
   )
{
   int i;

   SCIPdebugMessage("calling SCIPlpiChgBounds()\n");

   assert(lpi != NULL);
   assert(lpi->spx != NULL);
   assert(ind != NULL);
   assert(lb != NULL);
   assert(ub != NULL);

   invalidateSolution(lpi);

   assert( lpi->spx->preStrongbranchingBasisFreed() );

   try
   {
      for( i = 0; i < ncols; ++i )
      {
	 assert(0 <= ind[i] && ind[i] < lpi->spx->nCols());
	 lpi->spx->changeBounds(ind[i], lb[i], ub[i]);
         assert(lpi->spx->lower(ind[i]) <= lpi->spx->upper(ind[i]));
      }
   }
   catch(SPxException x)
   {
#ifndef NDEBUG
      std::string s = x.what();      
      SCIPmessagePrintWarning(lpi->messagehdlr, "SoPlex threw an exception: %s\n", s.c_str());
#endif
      return SCIP_LPERROR;
   }

   return SCIP_OKAY;
}

/** changes left and right hand sides of rows */
SCIP_RETCODE SCIPlpiChgSides(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   int                   nrows,              /**< number of rows to change sides for */
   const int*            ind,                /**< row indices */
   const SCIP_Real*      lhs,                /**< new values for left hand sides */
   const SCIP_Real*      rhs                 /**< new values for right hand sides */
   )
{
   int i;

   SCIPdebugMessage("calling SCIPlpiChgSides()\n");

   assert(lpi != NULL);
   assert(lpi->spx != NULL);
   assert(ind != NULL);
   assert(lhs != NULL);
   assert(rhs != NULL);

   invalidateSolution(lpi);

   assert( lpi->spx->preStrongbranchingBasisFreed() );

   try
   {
      for( i = 0; i < nrows; ++i )
      {
	 assert(0 <= ind[i] && ind[i] < lpi->spx->nRows());
	 lpi->spx->changeRange(ind[i], lhs[i], rhs[i]);
         assert(lpi->spx->lhs(ind[i]) <= lpi->spx->rhs(ind[i]));
      }
   }
   catch(SPxException x)
   {
#ifndef NDEBUG
      std::string s = x.what();      
      SCIPmessagePrintWarning(lpi->messagehdlr, "SoPlex threw an exception: %s\n", s.c_str());
#endif
      return SCIP_LPERROR;
   }

   return SCIP_OKAY;
}

/** changes a single coefficient */
SCIP_RETCODE SCIPlpiChgCoef(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   int                   row,                /**< row number of coefficient to change */
   int                   col,                /**< column number of coefficient to change */
   SCIP_Real             newval              /**< new value of coefficient */
   )
{
   SCIPdebugMessage("calling SCIPlpiChgCoef()\n");

   assert(lpi != NULL);
   assert(lpi->spx != NULL);
   assert(0 <= row && row < lpi->spx->nRows());
   assert(0 <= col && col < lpi->spx->nCols());

   invalidateSolution(lpi);

   assert( lpi->spx->preStrongbranchingBasisFreed() );

   SOPLEX_TRY( lpi->messagehdlr, lpi->spx->changeElement(row, col, newval) );

   return SCIP_OKAY;
}

/** changes the objective sense */
SCIP_RETCODE SCIPlpiChgObjsen(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   SCIP_OBJSEN           objsen              /**< new objective sense */
   )
{
   SCIPdebugMessage("calling SCIPlpiChgObjsen()\n");

   assert(lpi != NULL);
   assert(lpi->spx != NULL);

   invalidateSolution(lpi);

   assert( lpi->spx->preStrongbranchingBasisFreed() );

   SOPLEX_TRY( lpi->messagehdlr, lpi->spx->setSense(spxObjsen(objsen)) );

   return SCIP_OKAY;
}

/** changes objective values of columns in the LP */
SCIP_RETCODE SCIPlpiChgObj(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   int                   ncols,              /**< number of columns to change objective value for */
   int*                  ind,                /**< column indices to change objective value for */
   SCIP_Real*            obj                 /**< new objective values for columns */
   )
{
   int i;

   SCIPdebugMessage("calling SCIPlpiChgObj()\n");

   assert(lpi != NULL);
   assert(lpi->spx != NULL);
   assert(ind != NULL);
   assert(obj != NULL);

   invalidateSolution(lpi);

   assert( lpi->spx->preStrongbranchingBasisFreed() );

   try
   {
      for( i = 0; i < ncols; ++i )
      {
	 assert(0 <= ind[i] && ind[i] < lpi->spx->nCols());
	 lpi->spx->changeObj(ind[i], obj[i]);
      }
   }
   catch(SPxException x)
   {
#ifndef NDEBUG
      std::string s = x.what();      
      SCIPmessagePrintWarning(lpi->messagehdlr, "SoPlex threw an exception: %s\n", s.c_str());
#endif
      return SCIP_LPERROR;
   }

   return SCIP_OKAY;
}

/** multiplies a row with a non-zero scalar; for negative scalars, the row's sense is switched accordingly */
SCIP_RETCODE SCIPlpiScaleRow(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   int                   row,                /**< row number to scale */
   SCIP_Real             scaleval            /**< scaling multiplier */
   )
{
   SCIP_Real lhs;
   SCIP_Real rhs;

   SCIPdebugMessage("calling SCIPlpiScaleRow()\n");

   assert(lpi != NULL);
   assert(lpi->spx != NULL);
   assert(scaleval != 0.0);

   try 
   {
      invalidateSolution(lpi);

      assert( lpi->spx->preStrongbranchingBasisFreed() );

      /* get the row vector and the row's sides */
      SVector rowvec = lpi->spx->rowVector(row);
      lhs = lpi->spx->lhs(row);
      rhs = lpi->spx->rhs(row);

      /* scale the row vector */
      rowvec *= scaleval;

      /* adjust the sides */
      if( lhs > -soplex::infinity )
         lhs *= scaleval;
      else if( scaleval < 0.0 )
         lhs = soplex::infinity;
      if( rhs < soplex::infinity )
         rhs *= scaleval;
      else if( scaleval < 0.0 )
         rhs = -soplex::infinity;
      if( scaleval < 0.0 )
      {
         SCIP_Real oldlhs = lhs;
         lhs = rhs;
         rhs = oldlhs;
      }

      /* create the new row */
      LPRow lprow(lhs, rowvec, rhs);

      /* change the row in the LP */
      lpi->spx->changeRow(row, lprow);
      assert(lpi->spx->lhs(row) <= lpi->spx->rhs(row));
   }
   catch(SPxException x)
   {
#ifndef NDEBUG
      std::string s = x.what();      
      SCIPmessagePrintWarning(lpi->messagehdlr, "SoPlex threw an exception: %s\n", s.c_str());
#endif
      return SCIP_LPERROR;
   }

   return SCIP_OKAY;
}

/** multiplies a column with a non-zero scalar; the objective value is multiplied with the scalar, and the bounds
 *  are divided by the scalar; for negative scalars, the column's bounds are switched
 */
SCIP_RETCODE SCIPlpiScaleCol(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   int                   col,                /**< column number to scale */
   SCIP_Real             scaleval            /**< scaling multiplier */
   )
{
   SCIP_Real obj;
   SCIP_Real lb;
   SCIP_Real ub;

   SCIPdebugMessage("calling SCIPlpiScaleCol()\n");

   assert(lpi != NULL);
   assert(lpi->spx != NULL);
   assert(scaleval != 0.0);

   try
   {
      invalidateSolution(lpi);

      assert( lpi->spx->preStrongbranchingBasisFreed() );

      /* get the col vector and the col's bounds and objective value */
      SVector colvec = lpi->spx->colVector(col);
      obj = lpi->spx->obj(col);
      lb = lpi->spx->lower(col);
      ub = lpi->spx->upper(col);

      /* scale the col vector */
      colvec *= scaleval;

      /* scale the objective value */
      obj *= scaleval;

      /* adjust the bounds */
      if( lb > -soplex::infinity )
         lb /= scaleval;
      else if( scaleval < 0.0 )
         lb = soplex::infinity;
      if( ub < soplex::infinity )
         ub /= scaleval;
      else if( scaleval < 0.0 )
         ub = -soplex::infinity;
      if( scaleval < 0.0 )
      {
         SCIP_Real oldlb = lb;
         lb = ub;
         ub = oldlb;
      }

      /* create the new col (in LPCol's constructor, the upper bound is given first!) */
      LPCol lpcol(obj, colvec, ub, lb);

      /* change the col in the LP */
      lpi->spx->changeCol(col, lpcol);
      assert(lpi->spx->lower(col) <= lpi->spx->upper(col));
   }
   catch(SPxException x)
   {
#ifndef NDEBUG
      std::string s = x.what();      
      SCIPmessagePrintWarning(lpi->messagehdlr, "SoPlex threw an exception: %s\n", s.c_str());
#endif
      return SCIP_LPERROR;
   }

   return SCIP_OKAY;
}

/**@} */




/*
 * Data Accessing Methods
 */

/**@name Data Accessing Methods */
/**@{ */

/** gets the number of rows in the LP */
SCIP_RETCODE SCIPlpiGetNRows(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   int*                  nrows               /**< pointer to store the number of rows */
   )
{
   SCIPdebugMessage("calling SCIPlpiGetNRows()\n");

   assert(lpi != NULL);
   assert(lpi->spx != NULL);
   assert(nrows != NULL);

   *nrows = lpi->spx->nRows();

   return SCIP_OKAY;
}

/** gets the number of columns in the LP */
SCIP_RETCODE SCIPlpiGetNCols(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   int*                  ncols               /**< pointer to store the number of cols */
   )
{
   SCIPdebugMessage("calling SCIPlpiGetNCols()\n");

   assert(lpi != NULL);
   assert(lpi->spx != NULL);
   assert(ncols != NULL);

   *ncols = lpi->spx->nCols();

   return SCIP_OKAY;
}

/** gets the number of nonzero elements in the LP constraint matrix */
SCIP_RETCODE SCIPlpiGetNNonz(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   int*                  nnonz               /**< pointer to store the number of nonzeros */
   )
{
   int i;

   SCIPdebugMessage("calling SCIPlpiGetNNonz()\n");

   assert(lpi != NULL);
   assert(lpi->spx != NULL);
   assert(nnonz != NULL);

   /* SoPlex has no direct method to return the number of nonzeros, so we have to count them manually */
   *nnonz = 0;
   if( lpi->spx->nRows() < lpi->spx->nCols() )
   {
      for( i = 0; i < lpi->spx->nRows(); ++i )
         (*nnonz) += lpi->spx->rowVector(i).size();
   }
   else
   {
      for( i = 0; i < lpi->spx->nCols(); ++i )
         (*nnonz) += lpi->spx->colVector(i).size();
   }

   return SCIP_OKAY;
}

/** gets columns from LP problem object; the arrays have to be large enough to store all values
 *  Either both, lb and ub, have to be NULL, or both have to be non-NULL,
 *  either nnonz, beg, ind, and val have to be NULL, or all of them have to be non-NULL.
 */
SCIP_RETCODE SCIPlpiGetCols(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   int                   firstcol,           /**< first column to get from LP */
   int                   lastcol,            /**< last column to get from LP */
   SCIP_Real*            lb,                 /**< buffer to store the lower bound vector, or NULL */
   SCIP_Real*            ub,                 /**< buffer to store the upper bound vector, or NULL */
   int*                  nnonz,              /**< pointer to store the number of nonzero elements returned, or NULL */
   int*                  beg,                /**< buffer to store start index of each column in ind- and val-array, or NULL */
   int*                  ind,                /**< buffer to store column indices of constraint matrix entries, or NULL */
   SCIP_Real*            val                 /**< buffer to store values of constraint matrix entries, or NULL */
   )
{
   int i;
   int j;

   SCIPdebugMessage("calling SCIPlpiGetCols()\n");

   assert(lpi != NULL);
   assert(lpi->spx != NULL);
   assert(0 <= firstcol && firstcol <= lastcol && lastcol < lpi->spx->nCols());

   if( lb != NULL )
   {
      assert(ub != NULL);

      const Vector& lbvec = lpi->spx->lower();
      const Vector& ubvec = lpi->spx->upper();
      for( i = firstcol; i <= lastcol; ++i )
      {
         lb[i-firstcol] = lbvec[i];
         ub[i-firstcol] = ubvec[i];
      }
   }
   else
      assert(ub == NULL);

   if( nnonz != NULL )
   {
      *nnonz = 0;
      for( i = firstcol; i <= lastcol; ++i )
      {
         beg[i-firstcol] = *nnonz;
         const SVector& cvec = lpi->spx->colVector(i);
         for( j = 0; j < cvec.size(); ++j )
         {
            ind[*nnonz] = cvec.index(j);
            val[*nnonz] = cvec.value(j);
            (*nnonz)++;
         }
      }
   }
   else
   {
      assert(beg == NULL);
      assert(ind == NULL);
      assert(val == NULL);
   }

   return SCIP_OKAY;
}

/** gets rows from LP problem object; the arrays have to be large enough to store all values.
 *  Either both, lhs and rhs, have to be NULL, or both have to be non-NULL,
 *  either nnonz, beg, ind, and val have to be NULL, or all of them have to be non-NULL.
 */
SCIP_RETCODE SCIPlpiGetRows(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   int                   firstrow,           /**< first row to get from LP */
   int                   lastrow,            /**< last row to get from LP */
   SCIP_Real*            lhs,                /**< buffer to store left hand side vector, or NULL */
   SCIP_Real*            rhs,                /**< buffer to store right hand side vector, or NULL */
   int*                  nnonz,              /**< pointer to store the number of nonzero elements returned, or NULL */
   int*                  beg,                /**< buffer to store start index of each row in ind- and val-array, or NULL */
   int*                  ind,                /**< buffer to store row indices of constraint matrix entries, or NULL */
   SCIP_Real*            val                 /**< buffer to store values of constraint matrix entries, or NULL */
   )
{
   int i;
   int j;

   SCIPdebugMessage("calling SCIPlpiGetRows()\n");

   assert(lpi != NULL);
   assert(lpi->spx != NULL);
   assert(0 <= firstrow && firstrow <= lastrow && lastrow < lpi->spx->nRows());

   if( lhs != NULL )
   {
      assert(rhs != NULL);

      const Vector& lhsvec = lpi->spx->lhs();
      const Vector& rhsvec = lpi->spx->rhs();
      for( i = firstrow; i <= lastrow; ++i )
      {
         lhs[i-firstrow] = lhsvec[i];
         rhs[i-firstrow] = rhsvec[i];
      }
   }
   else
      assert(rhs == NULL);

   if( nnonz != NULL )
   {
      *nnonz = 0;
      for( i = firstrow; i <= lastrow; ++i )
      {
         beg[i-firstrow] = *nnonz;
         const SVector& rvec = lpi->spx->rowVector(i);
         for( j = 0; j < rvec.size(); ++j )
         {
            ind[*nnonz] = rvec.index(j);
            val[*nnonz] = rvec.value(j);
            (*nnonz)++;
         }
      }
   }
   else
   {
      assert(beg == NULL);
      assert(ind == NULL);
      assert(val == NULL);
   }

   return SCIP_OKAY;
}

/** gets column names */
SCIP_RETCODE SCIPlpiGetColNames(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   int                   firstcol,           /**< first column to get name from LP */
   int                   lastcol,            /**< last column to get name from LP */
   char**                colnames,           /**< pointers to column names (of size at least lastcol-firstcol+1) */
   char*                 namestorage,        /**< storage for col names */
   int                   namestoragesize,    /**< size of namestorage (if 0, storageleft returns the storage needed) */
   int*                  storageleft         /**< amount of storage left (if < 0 the namestorage was not big enough) */
   )
{
   assert( lpi != NULL );
   assert( lpi->spx != NULL );
   assert( colnames != NULL || namestoragesize == 0 );
   assert( namestorage != NULL || namestoragesize == 0 );
   assert( namestoragesize >= 0 );
   assert( storageleft != NULL );
   assert( 0 <= firstcol && firstcol <= lastcol && lastcol < lpi->spx->nCols() );

   SCIPdebugMessage("getting column names %d to %d\n", firstcol, lastcol);

   lpi->spx->getColNames(firstcol, lastcol, colnames, namestorage, namestoragesize, storageleft);

   return SCIP_OKAY;
}

/** gets row names */
SCIP_RETCODE SCIPlpiGetRowNames(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   int                   firstrow,           /**< first row to get name from LP */
   int                   lastrow,            /**< last row to get name from LP */
   char**                rownames,           /**< pointers to row names (of size at least lastrow-firstrow+1) */
   char*                 namestorage,        /**< storage for row names */
   int                   namestoragesize,    /**< size of namestorage (if 0, -storageleft returns the storage needed) */
   int*                  storageleft         /**< amount of storage left (if < 0 the namestorage was not big enough) */
   )
{
   assert( lpi != NULL );
   assert( lpi->spx != NULL );
   assert( rownames != NULL || namestoragesize == 0 );
   assert( namestorage != NULL || namestoragesize == 0 );
   assert( namestoragesize >= 0 );
   assert( storageleft != NULL );
   assert( 0 <= firstrow && firstrow <= lastrow && lastrow < lpi->spx->nRows() );

   SCIPdebugMessage("getting row names %d to %d\n", firstrow, lastrow);

   lpi->spx->getRowNames(firstrow, lastrow, rownames, namestorage, namestoragesize, storageleft);

   return SCIP_OKAY;
}

/** gets objective sense of the LP */
SCIP_RETCODE SCIPlpiGetObjsen(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   SCIP_OBJSEN*          objsen              /**< pointer to store objective sense */
   )
{
   SCIPdebugMessage("calling SCIPlpiGetObjsen()\n");

   assert(lpi != NULL);
   assert(lpi->spx != NULL);
   assert(objsen != NULL);

   *objsen = (lpi->spx->getSense() == SPxLP::MINIMIZE) ? SCIP_OBJSEN_MINIMIZE : SCIP_OBJSEN_MAXIMIZE;

   return SCIP_OKAY;
}

/** gets objective coefficients from LP problem object */
SCIP_RETCODE SCIPlpiGetObj(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   int                   firstcol,           /**< first column to get objective coefficient for */
   int                   lastcol,            /**< last column to get objective coefficient for */
   SCIP_Real*            vals                /**< array to store objective coefficients */
   )
{
   int i;

   SCIPdebugMessage("calling SCIPlpiGetObj()\n");

   assert(lpi != NULL);
   assert(lpi->spx != NULL);
   assert(0 <= firstcol && firstcol <= lastcol && lastcol < lpi->spx->nCols());
   assert(vals != NULL);
   
   for( i = firstcol; i <= lastcol; ++i )
      vals[i-firstcol] = lpi->spx->obj(i);

   return SCIP_OKAY;
}

/** gets current bounds from LP problem object */
SCIP_RETCODE SCIPlpiGetBounds(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   int                   firstcol,           /**< first column to get objective value for */
   int                   lastcol,            /**< last column to get objective value for */
   SCIP_Real*            lbs,                /**< array to store lower bound values, or NULL */
   SCIP_Real*            ubs                 /**< array to store upper bound values, or NULL */
   )
{
   int i;

   SCIPdebugMessage("calling SCIPlpiGetBounds()\n");

   assert(lpi != NULL);
   assert(lpi->spx != NULL);
   assert(0 <= firstcol && firstcol <= lastcol && lastcol < lpi->spx->nCols());
   
   for( i = firstcol; i <= lastcol; ++i )
   {
      if( lbs != NULL )
         lbs[i-firstcol] = lpi->spx->lower(i);
      if( ubs != NULL )
         ubs[i-firstcol] = lpi->spx->upper(i);
   }

   return SCIP_OKAY;
}

/** gets current row sides from LP problem object */
SCIP_RETCODE SCIPlpiGetSides(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   int                   firstrow,           /**< first row to get sides for */
   int                   lastrow,            /**< last row to get sides for */
   SCIP_Real*            lhss,               /**< array to store left hand side values, or NULL */
   SCIP_Real*            rhss                /**< array to store right hand side values, or NULL */
   )
{
   int i;

   SCIPdebugMessage("calling SCIPlpiGetSides()\n");

   assert(lpi != NULL);
   assert(lpi->spx != NULL);
   assert(0 <= firstrow && firstrow <= lastrow && lastrow < lpi->spx->nRows());
   
   for( i = firstrow; i <= lastrow; ++i )
   {
      if( lhss != NULL )
         lhss[i-firstrow] = lpi->spx->lhs(i);
      if( rhss != NULL )
         rhss[i-firstrow] = lpi->spx->rhs(i);
   }

   return SCIP_OKAY;
}

/** gets a single coefficient */
SCIP_RETCODE SCIPlpiGetCoef(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   int                   row,                /**< row number of coefficient */
   int                   col,                /**< column number of coefficient */
   SCIP_Real*            val                 /**< pointer to store the value of the coefficient */
   )
{
   SCIPdebugMessage("calling SCIPlpiGetCoef()\n");

   assert(lpi != NULL);
   assert(lpi->spx != NULL);
   assert(0 <= col && col < lpi->spx->nCols());
   assert(0 <= row && row < lpi->spx->nRows());
   assert(val != NULL);

   *val = lpi->spx->colVector(col)[row];

   return SCIP_OKAY;
}

/**@} */




/*
 * Solving Methods
 */

/**@name Solving Methods */
/**@{ */

/** solves LP -- used for both, primal and dual simplex, because SoPlex doesn't distinct the two cases */
static
SCIP_RETCODE spxSolve(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   SPxSolver::Representation rep,            /**< basis representation */
   SPxSolver::Type       type                /**< algorithm type */
   )
{
   assert( lpi != NULL );
   assert( lpi->spx != NULL );
   assert( rep == SPxSolver::ROW || rep == SPxSolver::COLUMN );
   assert( type == SPxSolver::ENTER || type == SPxSolver::LEAVE );

   SCIPdebugMessage("calling SoPlex solve(): %d cols, %d rows\n", lpi->spx->nCols(), lpi->spx->nRows());

   invalidateSolution(lpi);

   assert( lpi->spx->preStrongbranchingBasisFreed() );

   /* set basis representation and algorithm type */
   lpi->spx->setRep(rep);
   lpi->spx->setType(type);

#ifdef WITH_LPSCHECK
   lpi->spx->setDoubleCheck(CHECK_SPXSOLVE);
#endif

   SPxSolver::Status status = lpi->spx->solve();
   SCIPdebugMessage(" -> SoPlex status: %d, basis status: %d\n", lpi->spx->getStatus(), lpi->spx->basis().status());
   lpi->solved = TRUE;

   switch( status )
   {
   case SPxSolver::ABORT_TIME:
   case SPxSolver::ABORT_ITER:
   case SPxSolver::ABORT_VALUE:
   case SPxSolver::SINGULAR:
   case SPxSolver::REGULAR:
   case SPxSolver::UNKNOWN:
   case SPxSolver::OPTIMAL:
   case SPxSolver::UNBOUNDED:
   case SPxSolver::INFEASIBLE:
      return SCIP_OKAY;
   default:
      return SCIP_LPERROR;
   }  /*lint !e788*/
}

/** calls primal simplex to solve the LP */
SCIP_RETCODE SCIPlpiSolvePrimal(
   SCIP_LPI*             lpi                 /**< LP interface structure */
   )
{
   SCIPdebugMessage("calling SCIPlpiSolvePrimal()\n");

   SCIP_RETCODE retcode;
   SCIP_Bool rowrep;

   assert(lpi != NULL);
   assert(lpi->spx != NULL);

   /* first decide if we want to switch the basis representation; in order to avoid oscillatory behaviour, we add the
      factor 1.1 for switching back to column representation */
   if( lpi->rowrepswitch >= 0 )
   {
      rowrep = lpi->spx->rep() == SPxSolver::ROW;

      if( !rowrep )
         rowrep = lpi->spx->nRows() > lpi->spx->nCols() * (lpi->rowrepswitch);
      else
         rowrep = lpi->spx->nRows() * 1.1 > lpi->spx->nCols() * (lpi->rowrepswitch);
   }
   else
      rowrep = FALSE;

   /* SoPlex doesn't distinct between the primal and dual simplex; however
    * we can force SoPlex to start with the desired method:
    * If the representation is COLUMN:
    * - ENTER = PRIMAL 
    * - LEAVE = DUAL
    *
    * If the representation is ROW:
    * - ENTER = DUAL 
    * - LEAVE = PRIMAL
    */
   retcode = rowrep ? spxSolve(lpi, SPxSolver::ROW, SPxSolver::LEAVE) : spxSolve(lpi, SPxSolver::COLUMN, SPxSolver::ENTER);
   assert(!rowrep || lpi->spx->rep() == SPxSolver::ROW);
   assert(rowrep || lpi->spx->rep() == SPxSolver::COLUMN);

   return retcode;
}

/** calls dual simplex to solve the LP */
SCIP_RETCODE SCIPlpiSolveDual(
   SCIP_LPI*             lpi                 /**< LP interface structure */
   )
{
   SCIPdebugMessage("calling SCIPlpiSolveDual()\n");

   SCIP_RETCODE retcode;
   SCIP_Bool rowrep;

   assert(lpi != NULL);
   assert(lpi->spx != NULL);

   /* first decide if we want to switch the basis representation; in order to avoid oscillatory behaviour, we add the
      factor 1.1 for switching back to column representation */
   if( lpi->rowrepswitch >= 0 )
   {
      rowrep = lpi->spx->rep() == SPxSolver::ROW;

      if( !rowrep )
         rowrep = lpi->spx->nRows() > lpi->spx->nCols() * (lpi->rowrepswitch);
      else
         rowrep = lpi->spx->nRows() * 1.1 > lpi->spx->nCols() * (lpi->rowrepswitch);
   }
   else
      rowrep = FALSE;

   /* SoPlex doesn't distinct between the primal and dual simplex; however
    * we can force SoPlex to start with the desired method:
    * If the representation is COLUMN:
    * - ENTER = PRIMAL 
    * - LEAVE = DUAL
    *
    * If the representation is ROW:
    * - ENTER = DUAL 
    * - LEAVE = PRIMAL
    */
   retcode = rowrep ? spxSolve(lpi, SPxSolver::ROW, SPxSolver::ENTER) : spxSolve(lpi, SPxSolver::COLUMN, SPxSolver::LEAVE);
   assert(!rowrep || lpi->spx->rep() == SPxSolver::ROW);
   assert(rowrep || lpi->spx->rep() == SPxSolver::COLUMN);

   return retcode;
}

/** calls barrier or interior point algorithm to solve the LP with crossover to simplex basis */
SCIP_RETCODE SCIPlpiSolveBarrier(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   SCIP_Bool             crossover           /**< perform crossover */
   )
{  /*lint --e{715}*/
   SCIPdebugMessage("calling SCIPlpiSolveBarrier()\n");
   
   /* Since SoPlex does not support barrier we switch to DUAL */
   return SCIPlpiSolveDual(lpi);
}

/** start strong branching - call before any strongbranching */
SCIP_RETCODE SCIPlpiStartStrongbranch(
   SCIP_LPI*             lpi                 /**< LP interface structure */
   )
{
   assert( lpi->spx->preStrongbranchingBasisFreed() );
   lpi->spx->savePreStrongbranchingBasis();

   return SCIP_OKAY;
}

/** end strong branching - call after any strongbranching */
SCIP_RETCODE SCIPlpiEndStrongbranch(
   SCIP_LPI*             lpi                 /**< LP interface structure */
   )
{
   assert( ! lpi->spx->preStrongbranchingBasisFreed() );
   lpi->spx->restorePreStrongbranchingBasis();
   lpi->spx->freePreStrongbranchingBasis();

   return SCIP_OKAY;
}

/** performs strong branching iterations on one arbitrary candidate */
static
SCIP_RETCODE lpiStrongbranch(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   int                   col,                /**< column to apply strong branching on */
   SCIP_Real             psol,               /**< current primal solution value of column */
   int                   itlim,              /**< iteration limit for strong branchings */
   SCIP_Real*            down,               /**< stores dual bound after branching column down */
   SCIP_Real*            up,                 /**< stores dual bound after branching column up */
   SCIP_Bool*            downvalid,          /**< stores whether the returned down value is a valid dual bound;
                                              *   otherwise, it can only be used as an estimate value */
   SCIP_Bool*            upvalid,            /**< stores whether the returned up value is a valid dual bound;
                                              *   otherwise, it can only be used as an estimate value */
   int*                  iter                /**< stores total number of strong branching iterations, or -1; may be NULL */
   )
{
   SPxSCIP* spx;
   SPxSolver::Status status;
   SCIP_Real oldlb;
   SCIP_Real oldub;
   SCIP_Real newlb;
   SCIP_Real newub;
   bool fromparentbasis;
   bool error;
   int oldItlim;

   SCIPdebugMessage("calling SCIPlpiStrongbranch() on variable %d (%d iterations)\n", col, itlim);

   assert(lpi != NULL);
   assert(lpi->spx != NULL);
   assert(down != NULL);
   assert(up != NULL);
   assert(downvalid != NULL);
   assert(upvalid != NULL);

   spx = lpi->spx;
   status = SPxSolver::UNKNOWN;                      
   fromparentbasis = false;
   error = false;                                 
   oldItlim = spx->getIterationLimit();

   /* get current bounds of column */
   oldlb = spx->lower(col);
   oldub = spx->upper(col);

   *downvalid = FALSE;
   *upvalid = FALSE;

   if( iter != NULL )
      *iter = 0;

   /* set the algorithm type to use dual simplex */
   lpi->spx->setType( lpi->spx->rep() == SPxSolver::ROW ? SPxSolver::ENTER : SPxSolver::LEAVE);

   /* down branch */
   newub = EPSCEIL(psol-1.0, lpi->spx->feastol());
   if( newub >= oldlb - 0.5 )
   {
      SCIPdebugMessage("strong branching down on x%d (%g) with %d iterations\n", col, psol, itlim);

      spx->changeUpper(col, newub);
      assert(spx->lower(col) <= spx->upper(col));

      spx->setIterationLimit(itlim);
      do
      {
#ifdef WITH_LPSCHECK
         spx->setDoubleCheck(CHECK_SPXSTRONGBRANCH);
#endif
         status = spx->solve();
         SCIPdebugMessage(" --> Terminate with status %d\n", status);
         switch( status )
         {
         case SPxSolver::OPTIMAL:
            *down = spx->value();
            *downvalid = TRUE;
            SCIPdebugMessage(" --> Terminate with value %f\n", spx->value());
            break;
         case SPxSolver::ABORT_TIME: /* SoPlex does not return a proven dual bound, if it is aborted */
         case SPxSolver::ABORT_ITER:
         case SPxSolver::ABORT_CYCLING:
            *down = spx->value();
            break;
         case SPxSolver::ABORT_VALUE:
         case SPxSolver::INFEASIBLE:
            *down = spx->terminationValue();
            *downvalid = TRUE;
            break;
         default:
            error = true;
            break;
         }  /*lint !e788*/
         if( iter != NULL )
            (*iter) += spx->iterations();

#ifdef STRONGBRANCH_RESTOREBASIS
         /* we restore the pre-strong-branching basis by default (and don't solve again) */
         assert( ! spx->preStrongbranchingBasisFreed() );
         spx->restorePreStrongbranchingBasis();
         fromparentbasis = false;
#else
         /* if cycling or singular basis occured and we started not from the pre-strong-branching basis, then we restore the
          * pre-strong-branching basis and try again with reduced iteration limit */
         if( (status == SPxSolver::ABORT_CYCLING || status == SPxSolver::SINGULAR) && !fromparentbasis && spx->iterations() < itlim )
         {
            SCIPdebugMessage(" --> Repeat strong branching down with %d iterations after restoring basis\n", itlim - spx->iterations());
            spx->setIterationLimit(itlim - spx->iterations());
            assert( ! spx->hasPreStrongbranchingBasis() );
            spx->restorePreStrongbranchingBasis();
            fromparentbasis = true;
            error = false;
         }
         /* otherwise don't solve again */
         else
            fromparentbasis = false;
#endif
      }
      while( fromparentbasis );

      spx->changeUpper(col, oldub);
      assert(spx->lower(col) <= spx->upper(col));
   }
   else
   {
      *down = spx->terminationValue();
      *downvalid = TRUE;
   }

   /* up branch */
   if( !error )
   {
      newlb = EPSFLOOR(psol+1.0, lpi->spx->feastol());
      if( newlb <= oldub + 0.5 )
      {
         SCIPdebugMessage("strong branching  up  on x%d (%g) with %d iterations\n", col, psol, itlim);

         spx->changeLower(col, newlb);
         assert(spx->lower(col) <= spx->upper(col));

         spx->setIterationLimit(itlim);
         do
         {
#ifdef WITH_LPSCHECK
            spx->setDoubleCheck(CHECK_SPXSTRONGBRANCH);
#endif
            status = spx->solve();
            SCIPdebugMessage(" --> Terminate with status %d\n", status);
            switch( status )
            {
            case SPxSolver::OPTIMAL:
               *up = spx->value();
               *upvalid = TRUE;
               SCIPdebugMessage(" --> Terminate with value %f\n", spx->value());
               break;
            case SPxSolver::ABORT_TIME: /* SoPlex does not return a proven dual bound, if it is aborted */
            case SPxSolver::ABORT_ITER:
            case SPxSolver::ABORT_CYCLING:
               *up = spx->value();
               break;
            case SPxSolver::ABORT_VALUE:
            case SPxSolver::INFEASIBLE:
               *up = spx->terminationValue();
               *upvalid = TRUE;
               break;
            default:
               error = true;
               break;
            }  /*lint !e788*/
            if( iter != NULL )
               (*iter) += spx->iterations();

#ifdef STRONGBRANCH_RESTOREBASIS
            /* we restore the pre-strong-branching basis by default (and don't solve again) */
            assert( ! spx->preStrongbranchingBasisFreed() );
            spx->restorePreStrongbranchingBasis();
            fromparentbasis = false;
#else
            /* if cycling or singular basis occured and we started not from the pre-strong-branching basis, then we restore the
             * pre-strong-branching basis and try again with reduced iteration limit */
            else if( (status == SPxSolver::ABORT_CYCLING || status == SPxSolver::SINGULAR) && !fromparentbasis && spx->iterations() < itlim )
            {
               SCIPdebugMessage(" --> Repeat strong branching  up  with %d iterations after restoring basis\n", itlim - spx->iterations());
               assert( ! spx->hasPreStrongbranchingBasis() );
               spx->restorePreStrongbranchingBasis();
               spx->setIterationLimit(itlim - spx->iterations());
               error = false;
               fromparentbasis = true;
            }
            /* otherwise don't solve again */
            else
               fromparentbasis = false;
#endif
         }
         while( fromparentbasis );

         spx->changeLower(col, oldlb);
         assert(spx->lower(col) <= spx->upper(col));
      }
      else
      {
         *up = spx->terminationValue();
         *upvalid = TRUE;
      }
   }

   /* reset old iteration limit */
   spx->setIterationLimit(oldItlim);

   if( error )
   {
      SCIPdebugMessage("SCIPlpiStrongbranch() returned SoPlex status %d\n", int(status));
      return SCIP_LPERROR;
   }

   return SCIP_OKAY;
}

/** performs strong branching iterations on one @b fractional candidate */
SCIP_RETCODE SCIPlpiStrongbranchFrac(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   int                   col,                /**< column to apply strong branching on */
   SCIP_Real             psol,               /**< fractional current primal solution value of column */
   int                   itlim,              /**< iteration limit for strong branchings */
   SCIP_Real*            down,               /**< stores dual bound after branching column down */
   SCIP_Real*            up,                 /**< stores dual bound after branching column up */
   SCIP_Bool*            downvalid,          /**< stores whether the returned down value is a valid dual bound;
                                              *   otherwise, it can only be used as an estimate value */
   SCIP_Bool*            upvalid,            /**< stores whether the returned up value is a valid dual bound;
                                              *   otherwise, it can only be used as an estimate value */
   int*                  iter                /**< stores total number of strong branching iterations, or -1; may be NULL */
   )
{
   SCIP_RETCODE retcode;

   /* pass call on to lpiStrongbranch() */
   retcode = lpiStrongbranch(lpi, col, psol, itlim, down, up, downvalid, upvalid, iter);

   /* pass SCIP_LPERROR to SCIP without a back trace */
   if( retcode == SCIP_LPERROR )
      return SCIP_LPERROR;

   /* evaluate retcode */
   SCIP_CALL( retcode );

   return SCIP_OKAY;
}

/** performs strong branching iterations on given @b fractional candidates */
SCIP_RETCODE SCIPlpiStrongbranchesFrac(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   int*                  cols,               /**< columns to apply strong branching on */
   int                   ncols,              /**< number of columns */
   SCIP_Real*            psols,              /**< fractional current primal solution values of columns */
   int                   itlim,              /**< iteration limit for strong branchings */
   SCIP_Real*            down,               /**< stores dual bounds after branching columns down */
   SCIP_Real*            up,                 /**< stores dual bounds after branching columns up */
   SCIP_Bool*            downvalid,          /**< stores whether the returned down values are valid dual bounds;
                                              *   otherwise, they can only be used as an estimate values */
   SCIP_Bool*            upvalid,            /**< stores whether the returned up values are a valid dual bounds;
                                              *   otherwise, they can only be used as an estimate values */
   int*                  iter                /**< stores total number of strong branching iterations, or -1; may be NULL */
   )
{
   SCIP_RETCODE retcode;

   assert( iter != NULL );
   assert( cols != NULL );
   assert( psols != NULL );
   assert( down != NULL );
   assert( up != NULL );
   assert( downvalid != NULL );
   assert( upvalid != NULL );
   assert( down != NULL );

   if ( iter != NULL )
      *iter = 0;

   for (int j = 0; j < ncols; ++j)
   {
      /* pass call on to lpiStrongbranch() */
      retcode = lpiStrongbranch(lpi, cols[j], psols[j], itlim, &(down[j]), &(up[j]), &(downvalid[j]), &(upvalid[j]), iter);

      /* pass SCIP_LPERROR to SCIP without a back trace */
      if( retcode == SCIP_LPERROR )
         return SCIP_LPERROR;

      /* evaluate retcode */
      SCIP_CALL( retcode );
   }
   return SCIP_OKAY;
}

/** performs strong branching iterations on one candidate with @b integral value */
SCIP_RETCODE SCIPlpiStrongbranchInt(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   int                   col,                /**< column to apply strong branching on */
   SCIP_Real             psol,               /**< current integral primal solution value of column */
   int                   itlim,              /**< iteration limit for strong branchings */
   SCIP_Real*            down,               /**< stores dual bound after branching column down */
   SCIP_Real*            up,                 /**< stores dual bound after branching column up */
   SCIP_Bool*            downvalid,          /**< stores whether the returned down value is a valid dual bound;
                                              *   otherwise, it can only be used as an estimate value */
   SCIP_Bool*            upvalid,            /**< stores whether the returned up value is a valid dual bound;
                                              *   otherwise, it can only be used as an estimate value */
   int*                  iter                /**< stores total number of strong branching iterations, or -1; may be NULL */
   )
{
   SCIP_RETCODE retcode;

   /* pass call on to lpiStrongbranch() */
   retcode = lpiStrongbranch(lpi, col, psol, itlim, down, up, downvalid, upvalid, iter);

   /* pass SCIP_LPERROR to SCIP without a back trace */
   if( retcode == SCIP_LPERROR )
      return SCIP_LPERROR;

   /* evaluate retcode */
   SCIP_CALL( retcode );

   return SCIP_OKAY;
}

/** performs strong branching iterations on given candidates with @b integral values */
SCIP_RETCODE SCIPlpiStrongbranchesInt(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   int*                  cols,               /**< columns to apply strong branching on */
   int                   ncols,              /**< number of columns */
   SCIP_Real*            psols,              /**< current integral primal solution values of columns */
   int                   itlim,              /**< iteration limit for strong branchings */
   SCIP_Real*            down,               /**< stores dual bounds after branching columns down */
   SCIP_Real*            up,                 /**< stores dual bounds after branching columns up */
   SCIP_Bool*            downvalid,          /**< stores whether the returned down values are valid dual bounds;
                                              *   otherwise, they can only be used as an estimate values */
   SCIP_Bool*            upvalid,            /**< stores whether the returned up values are a valid dual bounds;
                                              *   otherwise, they can only be used as an estimate values */
   int*                  iter                /**< stores total number of strong branching iterations, or -1; may be NULL */
   )
{
   SCIP_RETCODE retcode;

   assert( iter != NULL );
   assert( cols != NULL );
   assert( psols != NULL );
   assert( down != NULL );
   assert( up != NULL );
   assert( downvalid != NULL );
   assert( upvalid != NULL );
   assert( down != NULL );

   if ( iter != NULL )
      *iter = 0;

   for (int j = 0; j < ncols; ++j)
   {
      /* pass call on to lpiStrongbranch() */
      retcode = lpiStrongbranch(lpi, cols[j], psols[j], itlim, &(down[j]), &(up[j]), &(downvalid[j]), &(upvalid[j]), iter);

      /* pass SCIP_LPERROR to SCIP without a back trace */
      if( retcode == SCIP_LPERROR )
         return SCIP_LPERROR;

      /* evaluate retcode */
      SCIP_CALL( retcode );
   }

   return SCIP_OKAY;
}
/**@} */




/*
 * Solution Information Methods
 */

/**@name Solution Information Methods */
/**@{ */

/** returns whether a solve method was called after the last modification of the LP */
SCIP_Bool SCIPlpiWasSolved(
   SCIP_LPI*             lpi                 /**< LP interface structure */
   )
{
   assert(lpi != NULL);

   return lpi->solved;
}

/** gets information about primal and dual feasibility of the current LP solution */
SCIP_RETCODE SCIPlpiGetSolFeasibility(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   SCIP_Bool*            primalfeasible,     /**< stores primal feasibility status */
   SCIP_Bool*            dualfeasible        /**< stores dual feasibility status */
   )
{
   SCIPdebugMessage("calling SCIPlpiGetSolFeasibility()\n");

   assert(lpi != NULL);
   assert(primalfeasible != NULL);
   assert(dualfeasible != NULL);

   *primalfeasible = SCIPlpiIsPrimalFeasible(lpi);
   *dualfeasible = SCIPlpiIsDualFeasible(lpi);

   return SCIP_OKAY;
}

/** returns TRUE iff LP is proven to have a primal unbounded ray (but not necessary a primal feasible point);
 *  this does not necessarily mean, that the solver knows and can return the primal ray
 */
SCIP_Bool SCIPlpiExistsPrimalRay(
   SCIP_LPI*             lpi                 /**< LP interface structure */
   )
{
   SCIPdebugMessage("calling SCIPlpiExistsPrimalRay()\n");

   assert(lpi != NULL);
   assert(lpi->spx != NULL);

   return (lpi->spx->getStatus() == SPxSolver::UNBOUNDED);
}

/** returns TRUE iff LP is proven to have a primal unbounded ray (but not necessary a primal feasible point),
 *  and the solver knows and can return the primal ray
 */
SCIP_Bool SCIPlpiHasPrimalRay(
   SCIP_LPI*             lpi                 /**< LP interface structure */
   )
{
   SCIPdebugMessage("calling SCIPlpiHasPrimalRay()\n");

   assert(lpi != NULL);
   assert(lpi->spx != NULL);

#if ((SOPLEX_VERSION == 150 && SOPLEX_SUBVERSION >= 2) || SOPLEX_VERSION > 150)
   return (lpi->spx->getStatus() == SPxSolver::UNBOUNDED);
#else
   return FALSE;
#endif
}

/** returns TRUE iff LP is proven to be primal unbounded */
SCIP_Bool SCIPlpiIsPrimalUnbounded(
   SCIP_LPI*             lpi                 /**< LP interface structure */
   )
{
   SCIPdebugMessage("calling SCIPlpiIsPrimalUnbounded()\n");

   assert(lpi != NULL);
   assert(lpi->spx != NULL);

   assert(lpi->spx->getStatus() != SPxSolver::UNBOUNDED || lpi->spx->basis().status() == SPxBasis::UNBOUNDED);

   /* if SoPlex returns unbounded, this may only mean that an unbounded ray is available, not necessarily a primal
    * feasible point; hence we have to check the perturbation
    */
   return (lpi->spx->getStatus() == SPxSolver::UNBOUNDED && !lpi->spx->isPerturbed());
}

/** returns TRUE iff LP is proven to be primal infeasible */
SCIP_Bool SCIPlpiIsPrimalInfeasible(
   SCIP_LPI*             lpi                 /**< LP interface structure */
   )
{
   SCIPdebugMessage("calling SCIPlpiIsPrimalInfeasible()\n");

   assert(lpi != NULL);
   assert(lpi->spx != NULL);

   return (lpi->spx->getStatus() == SPxSolver::INFEASIBLE);
}

/** returns TRUE iff LP is proven to be primal feasible */
SCIP_Bool SCIPlpiIsPrimalFeasible(
   SCIP_LPI*             lpi                 /**< LP interface structure */
   )
{
   SPxBasis::SPxStatus basestatus;

   SCIPdebugMessage("calling SCIPlpiIsPrimalFeasible()\n");

   assert(lpi != NULL);
   assert(lpi->spx != NULL);

   basestatus = lpi->spx->basis().status();

   /* note that the solver status may be ABORT_VALUE and the basis status optimal; if we are optimal, isPerturbed() may
    * still return true as long as perturbation plus violation is within tolerances
    */
   assert(basestatus == SPxBasis::OPTIMAL || lpi->spx->getStatus() != SPxSolver::OPTIMAL);

   return basestatus == SPxBasis::OPTIMAL ||
      ((basestatus == SPxBasis::PRIMAL || basestatus == SPxBasis::UNBOUNDED) && !lpi->spx->isPerturbed());
}

/** returns TRUE iff LP is proven to have a dual unbounded ray (but not necessary a dual feasible point);
 *  this does not necessarily mean, that the solver knows and can return the dual ray
 */
SCIP_Bool SCIPlpiExistsDualRay(
   SCIP_LPI*             lpi                 /**< LP interface structure */
   )
{
   SCIPdebugMessage("calling SCIPlpiExistsDualRay()\n");

   assert(lpi != NULL);
   assert(lpi->spx != NULL);

   return (lpi->spx->getStatus() == SPxSolver::INFEASIBLE);
}

/** returns TRUE iff LP is proven to have a dual unbounded ray (but not necessary a dual feasible point),
 *  and the solver knows and can return the dual ray
 */
SCIP_Bool SCIPlpiHasDualRay(
   SCIP_LPI*             lpi                 /**< LP interface structure */
   )
{
   SCIPdebugMessage("calling SCIPlpiHasDualRay()\n");

   assert(lpi != NULL);
   assert(lpi->spx != NULL);

   return (lpi->spx->getStatus() == SPxSolver::INFEASIBLE);
}

/** returns TRUE iff LP is dual unbounded */
SCIP_Bool SCIPlpiIsDualUnbounded(
   SCIP_LPI*             lpi                 /**< LP interface structure */
   )
{
   SCIPdebugMessage("calling SCIPlpiIsDualUnbounded()\n");

   assert(lpi != NULL);
   assert(lpi->spx != NULL);

   return (lpi->spx->getStatus() == SPxSolver::INFEASIBLE && lpi->spx->basis().status() == SPxBasis::DUAL
      && !lpi->spx->isPerturbed());
}

/** returns TRUE iff LP is dual infeasible */
SCIP_Bool SCIPlpiIsDualInfeasible(
   SCIP_LPI*             lpi                 /**< LP interface structure */
   )
{
   SCIPdebugMessage("calling SCIPlpiIsDualInfeasible()\n");

   assert(lpi != NULL);
   assert(lpi->spx != NULL);

   return (lpi->spx->getStatus() == SPxSolver::UNBOUNDED);
}

/** returns TRUE iff LP is proven to be dual feasible */
SCIP_Bool SCIPlpiIsDualFeasible(
   SCIP_LPI*             lpi                 /**< LP interface structure */
   )
{
   SCIPdebugMessage("calling SCIPlpiIsDualFeasible()\n");

   assert(lpi != NULL);
   assert(lpi->spx != NULL);

   /* note that the solver status may be ABORT_VALUE and the basis status optimal; if we are optimal, isPerturbed() may
    * still return true as long as perturbation plus violation is within tolerances
    */
   assert(lpi->spx->basis().status() == SPxBasis::OPTIMAL || lpi->spx->getStatus() != SPxSolver::OPTIMAL);

   return (lpi->spx->basis().status() == SPxBasis::OPTIMAL) ||
      (lpi->spx->basis().status() == SPxBasis::DUAL && !lpi->spx->isPerturbed());
}

/** returns TRUE iff LP was solved to optimality */
SCIP_Bool SCIPlpiIsOptimal(
   SCIP_LPI*             lpi                 /**< LP interface structure */
   )
{
   SCIPdebugMessage("calling SCIPlpiIsOptimal()\n");

   assert(lpi != NULL);
   assert((lpi->spx->basis().status() == SPxBasis::OPTIMAL)
      == (SCIPlpiIsPrimalFeasible(lpi) && SCIPlpiIsDualFeasible(lpi)));

   /* note that the solver status may be ABORT_VALUE and the basis status optimal; if we are optimal, isPerturbed() may
    * still return true as long as perturbation plus violation is within tolerances
    */
   return (lpi->spx->basis().status() == SPxBasis::OPTIMAL);
}

/** returns TRUE iff current LP basis is stable */
SCIP_Bool SCIPlpiIsStable(
   SCIP_LPI*             lpi                 /**< LP interface structure */
   )
{
   SCIPdebugMessage("calling SCIPlpiIsStable()\n");

   assert(lpi != NULL);
   assert(lpi->spx != NULL);

   return (lpi->spx->getStatus() != SPxSolver::ERROR && lpi->spx->getStatus() != SPxSolver::SINGULAR);
}

/** returns TRUE iff the objective limit was reached */
SCIP_Bool SCIPlpiIsObjlimExc(
   SCIP_LPI*             lpi                 /**< LP interface structure */
   )
{
   SCIPdebugMessage("calling SCIPlpiIsObjlimExc()\n");

   assert(lpi != NULL);
   assert(lpi->spx != NULL);

   return (lpi->spx->getStatus() == SPxSolver::ABORT_VALUE);
}

/** returns TRUE iff the iteration limit was reached */
SCIP_Bool SCIPlpiIsIterlimExc(
   SCIP_LPI*             lpi                 /**< LP interface structure */
   )
{
   SCIPdebugMessage("calling SCIPlpiIsIterlimExc()\n");

   assert(lpi != NULL);
   assert(lpi->spx != NULL);

   return (lpi->spx->getStatus() == SPxSolver::ABORT_ITER);
}

/** returns TRUE iff the time limit was reached */
SCIP_Bool SCIPlpiIsTimelimExc(
   SCIP_LPI*             lpi                 /**< LP interface structure */
   )
{
   SCIPdebugMessage("calling SCIPlpiIsTimelimExc()\n");

   assert(lpi != NULL);
   assert(lpi->spx != NULL);

   return (lpi->spx->getStatus() == SPxSolver::ABORT_TIME);
}

/** returns the internal solution status of the solver */
int SCIPlpiGetInternalStatus(
   SCIP_LPI*             lpi                 /**< LP interface structure */
   )
{
   SCIPdebugMessage("calling SCIPlpiIsTimelimExc()\n");

   assert(lpi != NULL);
   assert(lpi->spx != NULL);

   return static_cast<int>(lpi->spx->getStatus());
}

/** tries to reset the internal status of the LP solver in order to ignore an instability of the last solving call */
SCIP_RETCODE SCIPlpiIgnoreInstability(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   SCIP_Bool*            success             /**< pointer to store, whether the instability could be ignored */
   )
{  /*lint --e{715}*/
   SCIPdebugMessage("calling SCIPlpiIgnoreInstability()\n");

   assert(lpi != NULL);
   assert(lpi->spx != NULL);

   /* instable situations cannot be ignored */
   *success = FALSE;

   return SCIP_OKAY;
}

/** gets objective value of solution */
SCIP_RETCODE SCIPlpiGetObjval(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   SCIP_Real*            objval              /**< stores the objective value */
   )
{
   SCIPdebugMessage("calling SCIPlpiGetObjval()\n");

   assert(lpi != NULL);
   assert(lpi->spx != NULL);
   assert(objval != NULL);

   *objval = lpi->spx->value();

   return SCIP_OKAY;
}

/** gets primal and dual solution vectors */
SCIP_RETCODE SCIPlpiGetSol(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   SCIP_Real*            objval,             /**< stores the objective value, may be NULL if not needed */
   SCIP_Real*            primsol,            /**< primal solution vector, may be NULL if not needed */
   SCIP_Real*            dualsol,            /**< dual solution vector, may be NULL if not needed */
   SCIP_Real*            activity,           /**< row activity vector, may be NULL if not needed */
   SCIP_Real*            redcost             /**< reduced cost vector, may be NULL if not needed */
   )
{
   SCIPdebugMessage("calling SCIPlpiGetSol()\n");

   assert(lpi != NULL);
   assert(lpi->spx != NULL);

   if( objval != NULL )
      *objval = lpi->spx->value();

   try
   {
      if( primsol != NULL )
      {
         Vector tmp(lpi->spx->nCols(), primsol);
         (void)lpi->spx->getPrimal(tmp);
      }
      if( dualsol != NULL )
      {
         Vector tmp(lpi->spx->nRows(), dualsol);
         (void)lpi->spx->getDual(tmp);
      }
      if( activity != NULL )
      {
         Vector tmp(lpi->spx->nRows(), activity);
         (void)lpi->spx->getSlacks(tmp);  /* in SoPlex, the activities are called "slacks" */
      }
      if( redcost != NULL )
      {
         Vector tmp(lpi->spx->nCols(), redcost);
         (void)lpi->spx->getRedCost(tmp);
      }
   }
   catch(SPxException x)
   {
#ifndef NDEBUG
      std::string s = x.what();      
      SCIPmessagePrintWarning(lpi->messagehdlr, "SoPlex threw an exception: %s\n", s.c_str());
#endif
      return SCIP_LPERROR;
   }

   return SCIP_OKAY;
}

/** gets primal ray for unbounded LPs */
SCIP_RETCODE SCIPlpiGetPrimalRay(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   SCIP_Real*            ray                 /**< primal ray */
   )
{  /*lint --e{715}*/
   SCIPdebugMessage("calling SCIPlpiGetPrimalRay()\n");

   assert(lpi != NULL);
   assert(lpi->spx != NULL);

#if ((SOPLEX_VERSION == 150 && SOPLEX_SUBVERSION >= 2) || SOPLEX_VERSION > 150)
   try
   {
      Vector tmp(lpi->spx->nCols(), ray);
      (void)lpi->spx->getPrimalray(tmp);
   }
   catch(SPxException x)
   {
#ifndef NDEBUG
      std::string s = x.what();      
      SCIPmessagePrintWarning(lpi->messagehdlr, "SoPlex threw an exception: %s\n", s.c_str());
#endif
      return SCIP_LPERROR;
   }

   return SCIP_OKAY;
#else
   SCIPerrorMessage("SCIPlpiGetPrimalRay() not supported by SoPlex versions <= 1.5.0\n");
   return SCIP_LPERROR;
#endif
}

/** gets dual farkas proof for infeasibility */
SCIP_RETCODE SCIPlpiGetDualfarkas(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   SCIP_Real*            dualfarkas          /**< dual farkas row multipliers */
   )
{
   SCIPdebugMessage("calling SCIPlpiGetDualfarkas()\n");

   assert(lpi != NULL);
   assert(lpi->spx != NULL);

   try
   {
      Vector tmp(lpi->spx->nRows(), dualfarkas);
      (void)lpi->spx->getDualfarkas(tmp);
   }
   catch(SPxException x)
   {
#ifndef NDEBUG
      std::string s = x.what();      
      SCIPmessagePrintWarning(lpi->messagehdlr, "SoPlex threw an exception: %s\n", s.c_str());
#endif
      return SCIP_LPERROR;
   }

   return SCIP_OKAY;
}

/** gets the number of LP iterations of the last solve call */
SCIP_RETCODE SCIPlpiGetIterations(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   int*                  iterations          /**< pointer to store the number of iterations of the last solve call */
   )
{
   SCIPdebugMessage("calling SCIPlpiGetIterations()\n");

   assert(lpi != NULL);
   assert(lpi->spx != NULL);

   *iterations = lpi->spx->iterations();

   return SCIP_OKAY;
}

/** gets information about the quality of an LP solution
 *
 *  Such information is usually only available, if also a (maybe not optimal) solution is available.
 *  The LPI should return SCIP_INVALID for *quality, if the requested quantity is not available.
 */
SCIP_RETCODE SCIPlpiGetRealSolQuality(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   SCIP_LPSOLQUALITY     qualityindicator,   /**< indicates which quality should be returned */
   SCIP_Real*            quality             /**< pointer to store quality number */
   )
{
   SCIPdebugMessage("calling SCIPlpiGetRealSolQuality()\n");

   assert(lpi != NULL);
   assert(quality != NULL);

   *quality = SCIP_INVALID;

   return SCIP_OKAY;
}

/**@} */




/*
 * LP Basis Methods
 */

/**@name LP Basis Methods */
/**@{ */

/** Return reduced cost of column @c col if this is readily available, otherwise return 0.0 */
static
SCIP_RETCODE getRedCostEst(SPxSCIP* spx, int col, SCIP_Real* val)
{
   assert( spx != NULL );
   assert( val != NULL );

   *val = 0.0;

   /* Return if the vectors are not set up. The vectors are not set up if for instance we preformed
    * strong branching before. */
   if (! spx->isInitialized() )
      return SCIP_OKAY;

   assert( 0 <= col && col < spx->nCols() );

   if( spx->rep() == SPxSolver::COLUMN )
   {
      /* in column case the reduced costs are available: */
      if (spx->getSense() == SPxLP::MINIMIZE)
	 *val = spx->pVec()[col] - spx->maxObj()[col];
      else
	 *val = spx->maxObj()[col] - spx->pVec()[col];
   }
   else
   {
      assert( spx->rep() == SPxSolver::ROW );

      /* In row case for computing the reduced costs one needs to pass through the basis. We skip this expensive part. */
#if 0
      /* Here is the code necessary to compute the reduced costs for row representation: */
      SCIP_Real sign = 1.0;
      if ( spx->getSense() == SPxLP::MINIMIZE )
	 sign = -1.0;

      if ( spx->isColBasic(col) )
      {
	 /* It seems necessary to search through the basis in order to find the correct position */
         for (int i = spx->dim() - 1; i >= 0; --i)
         {
	    SPxId id = spx->basis().baseId(i);
            if ( id.isSPxColId() && col == spx->number(SPxColId(id)) )
	    {
	       *val = sign * spx->fVec()[i];
	       break;
	    }
	 }
      }
#endif
   }

   return SCIP_OKAY;
}



/** gets current basis status for columns and rows; arrays must be large enough to store the basis status */
SCIP_RETCODE SCIPlpiGetBase(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   int*                  cstat,              /**< array to store column basis status, or NULL */
   int*                  rstat               /**< array to store row basis status, or NULL */
   )
{
   int i;

   SCIPdebugMessage("calling SCIPlpiGetBase()\n");

   assert(lpi != NULL);
   assert(lpi->spx != NULL);

   assert( lpi->spx->preStrongbranchingBasisFreed() );

   if( rstat != NULL )
   {
      for( i = 0; i < lpi->spx->nRows(); ++i )
      {
         switch( lpi->spx->getBasisRowStatus(i) )
         {
         case SPxSolver::BASIC:
            rstat[i] = SCIP_BASESTAT_BASIC; /*lint !e641*/
            break;
         case SPxSolver::FIXED:
         case SPxSolver::ON_LOWER:
            rstat[i] = SCIP_BASESTAT_LOWER; /*lint !e641*/
            break;
         case SPxSolver::ON_UPPER:
            rstat[i] = SCIP_BASESTAT_UPPER; /*lint !e641*/
            break;
         case SPxSolver::ZERO:
            SCIPerrorMessage("slack variable has basis status ZERO (should not occur)\n");
            return SCIP_LPERROR;
         default:
            SCIPerrorMessage("invalid basis status\n");
            SCIPABORT();
            return SCIP_INVALIDDATA; /*lint !e527*/
         }
      }
   }

   if( cstat != NULL )
   {
      for( i = 0; i < lpi->spx->nCols(); ++i )
      {
	 SCIP_Real val = 0.0;
         switch( lpi->spx->getBasisColStatus(i) )
         {
         case SPxSolver::BASIC:
            cstat[i] = SCIP_BASESTAT_BASIC; /*lint !e641*/
            break;
         case SPxSolver::FIXED:
	    /* Get reduced cost estimation. If the estimation is not correct this should not hurt:
	     * If the basis is loaded into SoPlex again, the status is converted to FIXED again; in
	     * this case there is no problem at all. If the basis is saved and/or used in some other
	     * solver, it usually is very cheap to perform the pivots necessary to get an optimal
	     * basis. */
	    SCIP_CALL( getRedCostEst(lpi->spx, i, &val) );
	    if( val < 0.0 )  /* reduced costs < 0 => UPPER  else => LOWER */
	       cstat[i] = SCIP_BASESTAT_UPPER; /*lint !e641*/
	    else
	       cstat[i] = SCIP_BASESTAT_LOWER; /*lint !e641*/
            break;
         case SPxSolver::ON_LOWER:
            cstat[i] = SCIP_BASESTAT_LOWER; /*lint !e641*/
            break;
         case SPxSolver::ON_UPPER:
            cstat[i] = SCIP_BASESTAT_UPPER; /*lint !e641*/
            break;
         case SPxSolver::ZERO:
            cstat[i] = SCIP_BASESTAT_ZERO; /*lint !e641*/
            break;
         default:
            SCIPerrorMessage("invalid basis status\n");
            SCIPABORT();
            return SCIP_INVALIDDATA; /*lint !e527*/
         }
      }
   }

   return SCIP_OKAY;
}

/** sets current basis status for columns and rows */
SCIP_RETCODE SCIPlpiSetBase(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   int*                  cstat,              /**< array with column basis status */
   int*                  rstat               /**< array with row basis status */
   )
{
   int i;

   SCIPdebugMessage("calling SCIPlpiSetBase()\n");

   assert(lpi != NULL);
   assert(lpi->spx != NULL);
   assert(cstat != NULL || lpi->spx->nCols() == 0);
   assert(rstat != NULL || lpi->spx->nRows() == 0);

   assert( lpi->spx->preStrongbranchingBasisFreed() );
   invalidateSolution(lpi);

   SPxSolver::VarStatus* spxcstat = new SPxSolver::VarStatus[lpi->spx->nCols()];
   SPxSolver::VarStatus* spxrstat = new SPxSolver::VarStatus[lpi->spx->nRows()];

   for( i = 0; i < lpi->spx->nRows(); ++i )
   {
      switch( rstat[i] )
      {
      case SCIP_BASESTAT_LOWER:
         spxrstat[i] = SPxSolver::ON_LOWER;
         break;
      case SCIP_BASESTAT_BASIC:
         spxrstat[i] = SPxSolver::BASIC;
         break;
      case SCIP_BASESTAT_UPPER:
         spxrstat[i] = SPxSolver::ON_UPPER;
         break;
      case SCIP_BASESTAT_ZERO:
         SCIPerrorMessage("slack variable has basis status ZERO (should not occur)\n");
         delete[] spxcstat;
         delete[] spxrstat;
         return SCIP_LPERROR; /*lint !e429*/
      default:
         SCIPerrorMessage("invalid basis status\n");
         SCIPABORT();
         return SCIP_INVALIDDATA; /*lint !e527*/
      }
   }

   for( i = 0; i < lpi->spx->nCols(); ++i )
   {
      switch( cstat[i] )
      {
      case SCIP_BASESTAT_LOWER:
         spxcstat[i] = SPxSolver::ON_LOWER;
         break;
      case SCIP_BASESTAT_BASIC:
         spxcstat[i] = SPxSolver::BASIC;
         break;
      case SCIP_BASESTAT_UPPER:
         spxcstat[i] = SPxSolver::ON_UPPER;
         break;
      case SCIP_BASESTAT_ZERO:
         spxcstat[i] = SPxSolver::ZERO;
         break;
      default:
         SCIPerrorMessage("invalid basis status\n");
         SCIPABORT();
         return SCIP_INVALIDDATA; /*lint !e527*/
      }
   }

   SOPLEX_TRY( lpi->messagehdlr, lpi->spx->setBasis(spxrstat, spxcstat) );
   lpi->spx->updateStatus();

   delete[] spxcstat;
   delete[] spxrstat;

   return SCIP_OKAY;
}

/** returns the indices of the basic columns and rows */
SCIP_RETCODE SCIPlpiGetBasisInd(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   int*                  bind                /**< basic column n gives value n, basic row m gives value -1-m */
   )
{
   SCIPdebugMessage("calling SCIPlpiGetBasisInd()\n");

   assert(lpi != NULL);
   assert(lpi->spx != NULL);

   assert( lpi->spx->preStrongbranchingBasisFreed() );

   /* This function expects that the basis is a column basis. If SoPlex uses row representation, we
    * have to transform the basis and possibly compute a new factorization in the function
    * SCIPlpiGetBInvRow(), etc. */
   if( lpi->spx->rep() == SPxSolver::COLUMN )
   {
      /* for column representation, return the basis */
      SPxSolver* spx = lpi->spx;
      for (int i = 0; i < spx->nRows(); ++i)
      {
	 SPxId id = spx->basis().baseId(i);
	 if ( spx->isId(id) ) /* column id? */
	    bind[i] = spx->number(id);
	 else                 /* row id?    */
	    bind[i] = -1 - spx->number(id);
      }
   }
   else
   {
      SPxSolver* spx = lpi->spx;
      assert( spx->rep() == SPxSolver::ROW );

      /* for row representation, return the complement of the basis - need to loop through all rows and columns */
      int k = 0;
      for (int i = 0; i < spx->nRows(); ++i)
      {
	 if ( ! spx->isRowBasic(i) )
	    bind[k++] = -1 - i;
      }
      for (int j = 0; j < spx->nCols(); ++j)
      {
	 if ( ! spx->isColBasic(j) )
	    bind[k++] = j;
      }
      assert( k == spx->nRows() );
   }

   return SCIP_OKAY;
}


/* prepare a factorization of the basis matrix in column representation */
static
SCIP_RETCODE prepareFactorization(
   SCIP_LPI*   lpi
   )
{
   SCIPdebugMessage("Preparing factorization for computation of basis inverse.\n");

   try
   {
      /* if the factorization has not been set up, we compute a new factorization */
      if ( lpi->factorization == 0 )
      {
         SPxSolver* spx = lpi->spx;

         /* matrix to store columns */
         DataArray <const SVector*> matrix(spx->nRows());
      
         int k = 0;
         for (int i = 0; i < spx->nRows(); ++i)
         {
            if ( ! spx->isRowBasic(i) )
               matrix[k++] = new UnitVector(i);
         }
         for (int j = 0; j < spx->nCols(); ++j)
         {
            if ( ! spx->isColBasic(j) )
               matrix[k++] = &spx->colVector(j);
         }
         assert( k == spx->nRows() );
         assert( k == matrix.size() );

         /* compute factorization */
         lpi->factorization = new SLUFactor;
#ifndef NDEBUG
         SLinSolver::Status status = lpi->factorization->load(matrix.get_ptr(), k);
#else
         (void) lpi->factorization->load(matrix.get_ptr(), k);
#endif
         assert( status == SLinSolver::OK );
         assert( k == lpi->factorization->dim() );

         /* delete matrix columns corresponding to unit vectors */
         k = 0;
         for (int i = 0; i < spx->nRows(); ++i)
         {
            if ( ! spx->isRowBasic(i) )
               delete matrix[k++];
         }
      }
   }
   catch(SPxException x)
   {
#ifndef NDEBUG
      std::string s = x.what();      
      SCIPmessagePrintWarning(lpi->messagehdlr, "SoPlex threw an exception: %s\n", s.c_str());
#endif
      return SCIP_LPERROR;
   }
   return SCIP_OKAY;
}

/** get dense row of inverse basis matrix B^-1 */
SCIP_RETCODE SCIPlpiGetBInvRow(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   int                   r,                  /**< row number */
   SCIP_Real*            coef                /**< pointer to store the coefficients of the row */
   )
{
   SCIPdebugMessage("calling SCIPlpiGetBInvRow()\n");

   assert( lpi != NULL );
   assert( lpi->spx != NULL );
   assert( lpi->spx->preStrongbranchingBasisFreed() );

   try
   {
      SPxSolver* spx = lpi->spx;

      Vector x(spx->nRows(), coef); /* row of B^-1 has nrows entries - note that x is based on coef */
      DVector e(spx->nRows());      /* prepare unit vector */
      e.clear();
      e[r] = 1.0;

      /* in the column case use the existing factorization */
      if ( spx->rep() == SPxSolver::COLUMN )
      {
         assert( spx->dim() == spx->nRows() );
         assert( spx->coDim() == spx->nCols() );

         /* solve system "x = e_r^T * B^-1" to get r'th row of B^-1 */
         spx->basis().coSolve(x, e);
      }
      else
      {
         assert( spx->rep() == SPxSolver::ROW );
         assert( spx->dim() == spx->nCols() );
         assert( spx->coDim() == spx->nRows() );

         /* factorization is deleted in invalidateSolution() */
         SCIP_CALL( prepareFactorization(lpi) );
         assert( lpi->factorization != 0 );
         assert( lpi->factorization->dim() == spx->nRows() );
      
         /* solve system "x = e_r^T * B^-1" to get r'th row of B^-1 */
         lpi->factorization->solveLeft(x, e);
      }
   }
   catch(SPxException x)
   {
#ifndef NDEBUG
      std::string s = x.what();      
      SCIPmessagePrintWarning(lpi->messagehdlr, "SoPlex threw an exception: %s\n", s.c_str());
#endif
      return SCIP_LPERROR;
   }

   return SCIP_OKAY;
}

/** get dense column of inverse basis matrix B^-1 */
SCIP_RETCODE SCIPlpiGetBInvCol(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   int                   c,                  /**< column number of B^-1; this is NOT the number of the column in the LP;
                                              *   you have to call SCIPlpiGetBasisInd() to get the array which links the
                                              *   B^-1 column numbers to the row and column numbers of the LP!
                                              *   c must be between 0 and nrows-1, since the basis has the size
                                              *   nrows * nrows */
   SCIP_Real*            coef                /**< pointer to store the coefficients of the column */
   )
{
   SCIPdebugMessage("calling SCIPlpiGetBInvCol()\n");

   assert( lpi != NULL );
   assert( lpi->spx != NULL );
   assert( lpi->spx->preStrongbranchingBasisFreed() );

   try
   {
      SPxSolver* spx = lpi->spx;

      Vector x(spx->nRows(), coef); /* row of B^-1 has nrows entries - note that x is based on coef */
      DVector e(spx->nRows());      /* prepare unit vector */
      e.clear();
      e[c] = 1.0;

      /* in the column case use the existing factorization */
      if ( spx->rep() == SPxSolver::COLUMN )
      {
         assert( spx->dim() == spx->nRows() );
         assert( spx->coDim() == spx->nCols() );

         /* solve system "x = B^-1 * e_c" to get c'th column of B^-1 */
         spx->basis().solve(x, e);
      }
      else
      {
         assert( spx->rep() == SPxSolver::ROW );
         assert( spx->dim() == spx->nCols() );
         assert( spx->coDim() == spx->nRows() );
      
         /* factorization is deleted in invalidateSolution() */
         SCIP_CALL( prepareFactorization(lpi) );
         assert( lpi->factorization != 0 );
         assert( lpi->factorization->dim() == spx->nRows() );
      
         /* solve system "x = B^-1 * e_c" to get c'th column of B^-1 */
         lpi->factorization->solveRight(x, e);
      }
   }
   catch(SPxException x)
   {
#ifndef NDEBUG
      std::string s = x.what();      
      SCIPmessagePrintWarning(lpi->messagehdlr, "SoPlex threw an exception: %s\n", s.c_str());
#endif
      return SCIP_LPERROR;
   }

   return SCIP_OKAY;
}

/** get dense row of inverse basis matrix times constraint matrix B^-1 * A */
SCIP_RETCODE SCIPlpiGetBInvARow(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   int                   r,                  /**< row number */
   const SCIP_Real*      binvrow,            /**< row in (A_B)^-1 from prior call to SCIPlpiGetBInvRow(), or NULL */
   SCIP_Real*            coef                /**< vector to return coefficients */
   )
{
   SCIP_Real* buf;
   SCIP_Real* binv;
   int nrows;
   int ncols;
   int c;

   SCIPdebugMessage("calling SCIPlpiGetBInvARow()\n");

   assert(lpi != NULL);
   assert(lpi->spx != NULL);
   assert( lpi->spx->preStrongbranchingBasisFreed() );

   nrows = lpi->spx->nRows();
   ncols = lpi->spx->nCols();
   buf = NULL;

   /* get (or calculate) the row in B^-1 */
   if( binvrow == NULL )
   {
      SCIP_ALLOC( BMSallocMemoryArray(&buf, nrows) );
      SCIP_CALL( SCIPlpiGetBInvRow(lpi, r, buf) );
      binv = buf;
   }
   else
      binv = const_cast<SCIP_Real*>(binvrow);

   assert(binv != NULL);

   /* calculate the scalar product of the row in B^-1 and A */
   soplex::Vector binvvec(nrows, binv);
   for( c = 0; c < ncols; ++c )
      coef[c] = binvvec * lpi->spx->colVector(c);  /* scalar product */ /*lint !e1702*/

   /* free memory if it was temporarily allocated */
   BMSfreeMemoryArrayNull(&buf);

   return SCIP_OKAY;
}

/** get dense column of inverse basis matrix times constraint matrix B^-1 * A */
SCIP_RETCODE SCIPlpiGetBInvACol(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   int                   c,                  /**< column number */
   SCIP_Real*            coef                /**< vector to return coefficients */
   )
{
   SCIPdebugMessage("calling SCIPlpiGetBInvACol()\n");

   assert( lpi != NULL );
   assert( lpi->spx != NULL );
   assert( lpi->spx->preStrongbranchingBasisFreed() );

   try 
   {
      assert(lpi != NULL);
      assert(lpi->spx != NULL);
      SPxSolver* spx = lpi->spx;

      Vector x(spx->nRows(), coef); /* row of B^-1 has nrows entries - note that x is based on coef */
      DVector col(lpi->spx->nRows());

      /* extract column c of A */
      col = lpi->spx->colVector(c);

      /* in the column case use the existing factorization */
      if ( spx->rep() == SPxSolver::COLUMN )
      {
         assert( spx->dim() == spx->nRows() );
         assert( spx->coDim() == spx->nCols() );

         /* solve system "x = B^-1 * A_c" to get c'th column of B^-1 * A */
         lpi->spx->basis().solve(x, col);
      }
      else
      {
         assert( spx->rep() == SPxSolver::ROW );
         assert( spx->dim() == spx->nCols() );
         assert( spx->coDim() == spx->nRows() );

         /* factorization is deleted in invalidateSolution() */      
         SCIP_CALL( prepareFactorization(lpi) );
         assert( lpi->factorization != 0 );
         assert( lpi->factorization->dim() == spx->nRows() );

         /* solve system "x = B^-1 * A_c" to get c'th column of B^-1 * A */
         lpi->factorization->solveRight(x, col);
      }
   }
   catch(SPxException x)
   {
#ifndef NDEBUG
      std::string s = x.what();      
      SCIPmessagePrintWarning(lpi->messagehdlr, "SoPlex threw an exception: %s\n", s.c_str());
#endif
      return SCIP_LPERROR;
   }

   return SCIP_OKAY;
}

/**@} */




/*
 * LP State Methods
 */

/**@name LP State Methods */
/**@{ */

/** stores LPi state (like basis information) into lpistate object */
SCIP_RETCODE SCIPlpiGetState(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_LPISTATE**       lpistate            /**< pointer to LPi state information (like basis information) */
   )
{
   int ncols;
   int nrows;

   SCIPdebugMessage("calling SCIPlpiGetState()\n");

   assert(blkmem != NULL);
   assert(lpi != NULL);
   assert(lpi->spx != NULL);
   assert(lpistate != NULL);

   assert( lpi->spx->preStrongbranchingBasisFreed() );

   ncols = lpi->spx->nCols();
   nrows = lpi->spx->nRows();
   assert(ncols >= 0);
   assert(nrows >= 0);
   
   /* allocate lpistate data */
   SCIP_CALL( lpistateCreate(lpistate, blkmem, ncols, nrows) );

   /* allocate enough memory for storing uncompressed basis information */
   SCIP_CALL( ensureCstatMem(lpi, ncols) );
   SCIP_CALL( ensureRstatMem(lpi, nrows) );

   /* get unpacked basis information */
   SCIP_CALL( SCIPlpiGetBase(lpi, lpi->cstat, lpi->rstat) );

   /* pack LPi state data */
   (*lpistate)->ncols = ncols;
   (*lpistate)->nrows = nrows;
   lpistatePack(*lpistate, lpi->cstat, lpi->rstat);

   return SCIP_OKAY;
}

/** loads LPi state (like basis information) into solver; note that the LP might have been extended with additional
 *  columns and rows since the state was stored with SCIPlpiGetState()
 */
SCIP_RETCODE SCIPlpiSetState(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   BMS_BLKMEM*           /*blkmem*/,         /**< block memory */
   SCIP_LPISTATE*        lpistate            /**< LPi state information (like basis information) */
   )
{
   int lpncols;
   int lpnrows;
   int i;

   SCIPdebugMessage("calling SCIPlpiSetState()\n");

   assert(lpi != NULL);
   assert(lpi->spx != NULL);
   assert(lpistate != NULL);

   assert( lpi->spx->preStrongbranchingBasisFreed() );

   lpncols = lpi->spx->nCols();
   lpnrows = lpi->spx->nRows();
   assert(lpistate->ncols <= lpncols);
   assert(lpistate->nrows <= lpnrows);

   /* allocate enough memory for storing uncompressed basis information */
   SCIP_CALL( ensureCstatMem(lpi, lpncols) );
   SCIP_CALL( ensureRstatMem(lpi, lpnrows) );

   /* unpack LPi state data */
   lpistateUnpack(lpistate, lpi->cstat, lpi->rstat);

   /* extend the basis to the current LP beyond the previously existing columns */
   for( i = lpistate->ncols; i < lpncols; ++i )
   {
      SCIP_Real bnd = lpi->spx->lower(i);
      if ( SCIPlpiIsInfinity(lpi, REALABS(bnd)) )
      {
         /* if lower bound is +/- infinity -> try upper bound */
         bnd = lpi->spx->lower(i);
         if ( SCIPlpiIsInfinity(lpi, REALABS(bnd)) )
            lpi->cstat[i] = SCIP_BASESTAT_ZERO;  /* variable is free */
         else
            lpi->cstat[i] = SCIP_BASESTAT_UPPER; /* use finite upper bound */
      }
      else
         lpi->cstat[i] = SCIP_BASESTAT_LOWER;    /* use finite lower bound */
   }
   for( i = lpistate->nrows; i < lpnrows; ++i )
      lpi->rstat[i] = SCIP_BASESTAT_BASIC; /*lint !e641*/

   /* load basis information */
   SCIP_CALL( SCIPlpiSetBase(lpi, lpi->cstat, lpi->rstat) );

   return SCIP_OKAY;
}

/** clears current LPi state (like basis information) of the solver */
SCIP_RETCODE SCIPlpiClearState(
   SCIP_LPI*             lpi                 /**< LP interface structure */
   )
{  /*lint --e{715}*/
   SCIPdebugMessage("calling SCIPlpiClearState()\n");

   assert(lpi != NULL);
   assert(lpi->spx != NULL);

   try
   {
      lpi->spx->reLoad();
   }
   catch(SPxException x)
   {
#ifndef NDEBUG
      std::string s = x.what();
      SCIPmessagePrintWarning(lpi->messagehdlr, "SoPlex threw an exception: %s\n", s.c_str());
#endif
      assert( lpi->spx->getStatus() != SPxSolver::OPTIMAL );
      return SCIP_LPERROR;
   }

   return SCIP_OKAY;
}

/** frees LPi state information */
SCIP_RETCODE SCIPlpiFreeState(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_LPISTATE**       lpistate            /**< pointer to LPi state information (like basis information) */
   )
{  /*lint --e{715}*/
   SCIPdebugMessage("calling SCIPlpiFreeState()\n");

   assert(lpi != NULL);
   assert(lpistate != NULL);

   if ( *lpistate != NULL )
      lpistateFree(lpistate, blkmem);

   return SCIP_OKAY;
}

/** checks, whether the given LP state contains simplex basis information */
SCIP_Bool SCIPlpiHasStateBasis(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   SCIP_LPISTATE*        lpistate            /**< LP state information (like basis information) */
   )
{  /*lint --e{715}*/
   return TRUE;
}

/** reads LP state (like basis information from a file */
SCIP_RETCODE SCIPlpiReadState(
   SCIP_LPI*             lpi,               /**< LP interface structure */
   const char*           fname              /**< file name */
   )
{
   SCIPdebugMessage("calling SCIPlpiReadState()\n");

   assert( lpi->spx->preStrongbranchingBasisFreed() );

   bool success;
   SOPLEX_TRY( lpi->messagehdlr, success = lpi->spx->readBasisFile(fname, 0, 0) );

   return success ? SCIP_OKAY : SCIP_LPERROR;
}

/** writes LP state (like basis information) to a file */
SCIP_RETCODE SCIPlpiWriteState(
   SCIP_LPI*             lpi,            /**< LP interface structure */
   const char*           fname           /**< file name */
   )
{
   SCIPdebugMessage("calling SCIPlpiWriteState()\n");

   assert( lpi->spx->preStrongbranchingBasisFreed() );

   bool res;
   SOPLEX_TRY( lpi->messagehdlr, res = lpi->spx->writeBasisFile(fname, 0, 0) );

   if ( ! res )
      return SCIP_LPERROR;

   return SCIP_OKAY;
}

/**@} */




/*
 * Parameter Methods
 */

/**@name Parameter Methods */
/**@{ */

/** gets integer parameter of LP */
SCIP_RETCODE SCIPlpiGetIntpar(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   SCIP_LPPARAM          type,               /**< parameter number */
   int*                  ival                /**< buffer to store the parameter value */
   )
{
   SCIPdebugMessage("calling SCIPlpiGetIntpar()\n");

   assert(lpi != NULL);
   assert(lpi->spx != NULL);
   assert(ival != NULL);

   switch( type )
   {
   case SCIP_LPPAR_FROMSCRATCH:
      *ival = lpi->spx->getFromScratch();
      break;
   case SCIP_LPPAR_LPINFO:
      *ival = lpi->spx->getLpInfo();
      break;
   case SCIP_LPPAR_LPITLIM:
      *ival = lpi->spx->getIterationLimit();
      break;
   case SCIP_LPPAR_PRESOLVING:
      *ival = lpi->spx->getPresolving();
      break;
   case SCIP_LPPAR_PRICING:
      *ival = (int)lpi->pricing;
      break;
   case SCIP_LPPAR_SCALING:
      *ival = lpi->spx->getScaling();
      break;
   default:
      return SCIP_PARAMETERUNKNOWN;
   }  /*lint !e788*/

   return SCIP_OKAY;
}

/** sets integer parameter of LP */
SCIP_RETCODE SCIPlpiSetIntpar(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   SCIP_LPPARAM          type,               /**< parameter number */
   int                   ival                /**< parameter value */
   )
{
   SCIPdebugMessage("calling SCIPlpiSetIntpar()\n");

   assert(lpi != NULL);
   assert(lpi->spx != NULL);

   switch( type )
   {
   case SCIP_LPPAR_FROMSCRATCH:
      assert(ival == TRUE || ival == FALSE);
      lpi->spx->setFromScratch(bool(ival));
      break;
   case SCIP_LPPAR_LPINFO:
      assert(ival == TRUE || ival == FALSE);
      lpi->spx->setLpInfo(bool(ival));
      break;
   case SCIP_LPPAR_LPITLIM:
      assert(ival >= -1);
      lpi->spx->setIterationLimit(ival);
      break;
   case SCIP_LPPAR_PRESOLVING:
      assert(ival == TRUE || ival == FALSE);
      lpi->spx->setPresolving(bool(ival));
      break;
   case SCIP_LPPAR_PRICING:
      lpi->pricing = (SCIP_PRICING)ival;
      switch( lpi->pricing )
      {
      case SCIP_PRICING_LPIDEFAULT:
      case SCIP_PRICING_AUTO:
         lpi->spx->setAutoPricer();
         break;
      case SCIP_PRICING_FULL:
         lpi->spx->setFullPricer();
         break;
      case SCIP_PRICING_PARTIAL:
         lpi->spx->setParmultPricer();
         break;
      case SCIP_PRICING_STEEP:
         lpi->spx->setSteepPricer();
	 break;
      case SCIP_PRICING_STEEPQSTART:
         lpi->spx->setSteepQStartPricer();
	 break;
      case SCIP_PRICING_DEVEX:
         lpi->spx->setDevexPricer();
	 break;
      default:
         return SCIP_LPERROR;
      }
      break;
   case SCIP_LPPAR_SCALING:
      assert(ival == TRUE || ival == FALSE);
      lpi->spx->setScaling(bool(ival));
      break;
   default:
      return SCIP_PARAMETERUNKNOWN;
   }  /*lint !e788*/

   return SCIP_OKAY;
}

/** gets floating point parameter of LP */
SCIP_RETCODE SCIPlpiGetRealpar(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   SCIP_LPPARAM          type,               /**< parameter number */
   SCIP_Real*            dval                /**< buffer to store the parameter value */
   )
{
   SCIPdebugMessage("calling SCIPlpiGetRealpar()\n");

   assert(lpi != NULL);
   assert(lpi->spx != NULL);
   assert(dval != NULL);

   switch( type )
   {
   case SCIP_LPPAR_FEASTOL:
      *dval = lpi->spx->feastol();
      break;
#if ((SOPLEX_VERSION == 160 && SOPLEX_SUBVERSION >= 5) || SOPLEX_VERSION > 160)
   case SCIP_LPPAR_DUALFEASTOL:
      *dval = lpi->spx->opttol();
      break;
#endif
   case SCIP_LPPAR_LOBJLIM:
      *dval = lpi->spx->getObjLoLimit();
      break;
   case SCIP_LPPAR_UOBJLIM:
      *dval = lpi->spx->getObjUpLimit();
      break;
   case SCIP_LPPAR_LPTILIM:
      *dval = lpi->spx->terminationTime();
      break;
   case SCIP_LPPAR_ROWREPSWITCH:
      *dval = lpi->rowrepswitch;
      break;
   default:
      return SCIP_PARAMETERUNKNOWN;
   }  /*lint !e788*/

   return SCIP_OKAY;
}

/** sets floating point parameter of LP */
SCIP_RETCODE SCIPlpiSetRealpar(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   SCIP_LPPARAM          type,               /**< parameter number */
   SCIP_Real             dval                /**< parameter value */
   )
{
   SCIPdebugMessage("calling SCIPlpiSetRealpar()\n");

   assert(lpi != NULL);
   assert(lpi->spx != NULL);

   switch( type )
   {
   case SCIP_LPPAR_FEASTOL:
      lpi->spx->setFeastol(dval);
      break;
#if ((SOPLEX_VERSION == 160 && SOPLEX_SUBVERSION >= 5) || SOPLEX_VERSION > 160)
   case SCIP_LPPAR_DUALFEASTOL:
      lpi->spx->setOpttol(dval);
      break;
#endif
   case SCIP_LPPAR_LOBJLIM:
      lpi->spx->setObjLoLimit(dval);
      break;
   case SCIP_LPPAR_UOBJLIM:
      lpi->spx->setObjUpLimit(dval);
      break;
   case SCIP_LPPAR_LPTILIM:
      lpi->spx->setTerminationTime(dval);
      break;
   case SCIP_LPPAR_ROWREPSWITCH:
      assert(dval >= -1.5);
      lpi->rowrepswitch = dval;
      break;
   default:
      return SCIP_PARAMETERUNKNOWN;
   }  /*lint !e788*/

   return SCIP_OKAY;
}

/**@} */




/*
 * Numerical Methods
 */

/**@name Numerical Methods */
/**@{ */

/** returns value treated as infinity in the LP solver */
SCIP_Real SCIPlpiInfinity(
   SCIP_LPI*             /*lpi*/             /**< LP interface structure */
   )
{
   SCIPdebugMessage("calling SCIPlpiInfinity()\n");

   return soplex::infinity;
}

/** checks if given value is treated as infinity in the LP solver */
SCIP_Bool SCIPlpiIsInfinity(
   SCIP_LPI*             /*lpi*/,            /**< LP interface structure */
   SCIP_Real             val
   )
{
   SCIPdebugMessage("calling SCIPlpiIsInfinity()\n");

   return (val >= soplex::infinity);
}

/**@} */




/*
 * File Interface Methods
 */

/**@name File Interface Methods */
/**@{ */

/** returns, whether the given file exists */
static
SCIP_Bool fileExists(
   const char*           filename            /**< file name */
   )
{
   FILE* f;

   f = fopen(filename, "r");
   if( f == NULL )
      return FALSE;

   fclose(f);

   return TRUE;
}

/** reads LP from a file */
SCIP_RETCODE SCIPlpiReadLP(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   const char*           fname               /**< file name */
   )
{
   SCIPdebugMessage("calling SCIPlpiReadLP()\n");

   assert(lpi != NULL);
   assert(lpi->spx != NULL);

   assert( lpi->spx->preStrongbranchingBasisFreed() );

   if( !fileExists(fname) )
      return SCIP_NOFILE;

   try
   {
      if( !lpi->spx->readLP(fname) )
	 return SCIP_READERROR;
   }
   catch(SPxException x)
   {
#ifndef NDEBUG
      std::string s = x.what();      
      SCIPmessagePrintWarning(lpi->messagehdlr, "SoPlex threw an exception: %s\n", s.c_str());
#endif
      return SCIP_READERROR;
   }
   
   return SCIP_OKAY;
}

/** writes LP to a file */
SCIP_RETCODE SCIPlpiWriteLP(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   const char*           fname               /**< file name */
   )
{
   SCIPdebugMessage("calling SCIPlpiWriteLP()\n");

   assert(lpi != NULL);
   assert(lpi->spx != NULL);

   try
   {
      lpi->spx->writeFile(fname);
   }
   catch(SPxException x)
   {
#ifndef NDEBUG
      std::string s = x.what();
      SCIPmessagePrintWarning(lpi->messagehdlr, "SoPlex threw an exception: %s\n", s.c_str());
#endif
      return SCIP_WRITEERROR;
   }

   return SCIP_OKAY;
}

/**@} */

