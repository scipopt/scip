/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2003 Tobias Achterberg                              */
/*                            Thorsten Koch                                  */
/*                  2002-2003 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the SCIP Academic Licence.        */
/*                                                                           */
/*  You should have received a copy of the SCIP Academic License             */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   lp_spx.cpp
 * @brief  LP interface for SOPLEX
 * @author Tobias Achterberg
 */

/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/


#if 0 /* ??? TODO */


#include <fstream>
#include <stdio.h>
#include <string.h>
#include <assert.h>

#include "spxdefines.h"
#include "soplex.h"
#include "spxlp.h"
#include "slufactor.h"
#include "spxsteeppr.h"
#include "spxfastrt.h"
#include "nameset.h"
#include "didxset.h"

/* set DEBUGGING to 3 if you want to be informed of every LP interface call */
#define DEBUGGING 0
/* set DEBUG_FILEOUTPUT if you want to get LP files dumped */
/*#define DEBUG_FILEOUTPUT  1*/


extern "C" 
{
#include "func.h"
#include "lp.h"
#include "def.h"
#include "struct.h"
#include "misc.h"
}

/* SOPLEX doesn't have an environment.
 * All information is stored in the SOPLEX object.
 */
static int dummy_environment;

using namespace soplex;

#if DEBUGGING >= 2
extern "C" void print_base(
   char dir, int cols, int rows, int* cstat, int* rstat)
{
   static char ctext[] = { 'L', 'B', 'U', 'F' };
   static char rtext[] = { 'l', 'b', 'u', 'f' };

   int i;

   std::cerr << "\nC " << dir << " ";
   if (cstat != 0)
      for(i = 0; i < cols; i++)
         std::cerr << ctext[cstat[i]];
   std::cerr << "\nR " << dir << " ";
   if (rstat != 0)
      for(i = 0; i < rows; i++)
         std::cerr << rtext[rstat[i]];
   std::cerr << "\n";
}
#endif /*DEBUGGING >= 2*/

class SPxSIP : public SoPlex
{
   SLUFactor    m_slu;
   SPxSteepPR   m_price;
   SPxFastRT    m_ratio;
   char*        m_probname;
   bool         m_fromscratch;  /**< use old basis indicator */
   double       m_objLoLimit;
   double       m_objUpLimit;
   Status       m_stat;
   NameSet      m_colnames;
   NameSet      m_rownames;

public:
   static SPxSIP* fromSIPLP(SIPLP lp);

   bool getFromScratch() const
   {
      return m_fromscratch;
   }
   /**@todo If m_fromscratch is true, a new (slack/starter) basis
    *       should be generated.
    */
   void setFromScratch(bool fs)
   {
      m_fromscratch = fs;
   }

   /* Was macht das? Brauchen wir das?
    * void splitLP()
    * {
    *   subcovectors.reSize( 1 ) ;
    * }
    */

   SPxSIP() 
      : SoPlex(LEAVE, COLUMN)
      , m_probname(0)
      , m_fromscratch(false)
      , m_objLoLimit(-SIP_INFINITY)
      , m_objUpLimit(SIP_INFINITY)
      , m_stat(NO_PROBLEM)
   {
      setSolver(&m_slu);
      setTester(&m_ratio);
      setPricer(&m_price);
      /* no starter, no simplifier */

      Param::setVerbose(0);
      m_slu.setUtype(SLUFactor::ETA);
   }
   virtual ~SPxSIP()
   {
      if (m_probname != NULL)
         spx_free(m_probname);  /*lint !e1551*/
   }
   void setProbname(const char* probname)
   {
      assert(probname != NULL);
      if (m_probname != NULL)
         spx_free(m_probname);
      spx_alloc(m_probname, strlen(probname) + 1);
      strcpy(m_probname, probname);
   }
   void setObjLoLimit(double limit)
   {
      m_objLoLimit = limit;
   }
   void setObjUpLimit(double limit)
   {
      m_objUpLimit = limit;
   }
   virtual Status solve()
   {
      if( getFromScratch() )
      {
         SoPlex::reLoad();
         setFromScratch( false );
      }
      m_stat = SoPlex::solve();

      assert(rep() == COLUMN);

      if (m_stat == OPTIMAL)
      {
         double objval = value();

         if ((objval > m_objUpLimit) || (objval < m_objLoLimit))
            m_stat = ABORT_VALUE;
      }
      return m_stat;
   }
   Status getStatus() const
   {
      return m_stat;
   }

   /* Problem, es ist nicht moeglich bei einem NameSet sowohl den
    * Namen als auch den DataKey zu uebergeben. Das waere aber die
    * korrekte Implementierung zumal die LPxLP::addRows usw. 
    * jeweils auch eine Variante haben, die die neuen ID's 
    * zurueckliefert.
    * Dafuer brauchen wir aber erst eine neue NameSet Implementierung.
    * Das hier unten klappt nicht, weil nach removeRows aus irgendweinem
    * dummen Grund die ID's fuer das NameSet anders vergeben werden als 
    * fuer das RowSet.
    */
   void addRowsWithNames(const LPRowSet& rset, const NameSet& rnames)
   {
      assert(m_rownames.num() == nRows());
      assert(rnames.num() == rset.num());

#if 0
      std::cerr << "[AR " << m_rownames.num() << " " 
                << nRows() << " + "
                << rnames.num() << std::endl;

      int i; 

      for(i = 0; i < rset.num(); i++)
      {
         if (m_rownames.has(rnames[i]))
         {
            std::cerr << "*** Doublicate " << rnames[i] << std::endl
                      << m_rownames;
            abort();
         }
      }
      SPxRowId* newIdsR = new SPxRowId[rset.num()];
      SPxRowId* newIdsN = new SPxRowId[rset.num()];
      addRows(newIdsR, rset);
      m_rownames.add(newIdsN, rnames);

      for(i = 0; i < rset.num(); i++)
         std::cerr << "AR " << newIdsR[i].getIdx() << " "
                   << newIdsN[i].getIdx() << " "
                   << rnames[i] << " / "
                   << m_rownames[newIdsN[i]] 
                   << std::endl;

      for(i = 0; i < rset.num(); i++)
         assert(newIdsR[i].getIdx() ==  newIdsN[i].getIdx());

      delete [] newIdsR;
      delete [] newIdsN;

      std::cerr << "]AR " << m_rownames.num() << " " 
                << nRows() << " + "
                << rnames.num() << std::endl;
#endif

#ifndef NDEBUG      
      for(int i = 0; i < rset.num(); i++)
         if (m_rownames.has(rnames[i]))
            std::cerr << "*** Doublicate rowname " << rnames[i] << std::endl << m_rownames;
#endif /*NDEBUG*/
      addRows(rset);
      m_rownames.add(rnames);

      assert(m_rownames.num() == nRows());
   }
   void addColsWithNames(const LPColSet& cset, const NameSet& cnames)
   {
      assert(m_colnames.num() == nCols());
      assert(cnames.num() == cset.num());

      addCols(cset);
      m_colnames.add(cnames);

      assert(m_colnames.num() == nCols());
   }
   virtual void removeRowsWithNames(int dstat[])
   {
#if 0
      int count = 0;
#endif
      assert(m_rownames.num() == nRows());

      /* SOPLEX marks rows to be deleted with a value < 0,
       * while CPLEX marks them with dstat[i] == 1 
       */
      for(int i = 0; i < nRows(); i++)
      {
         if (dstat[i] != 0)
         {
            dstat[i] = -1;
#if 0
            count++;
            std::cerr << "-RR " << rId(i).getIdx() << std::endl;
#endif
         }
      }
      removeRows( dstat );
      m_rownames.remove( dstat );
#if 0
      std::cerr << m_rownames;

      std::cerr << "*RR " << m_rownames.num() << " " 
                << nRows() << " - "
                << count << std::endl;
#endif
      assert(m_rownames.num() == nRows());
   }
   virtual void clearWithNames()
   {
      clear();
      m_rownames.clear();
      m_colnames.clear();

      if (m_probname != NULL)
         spx_free(m_probname);  /*lint !e1551*/

      m_stat = NO_PROBLEM;
   }
   void dumpBasis(const char* filename)
   {
      std::ofstream file(filename);

      if (file.good())
         writeBasis(file, m_rownames, m_colnames);
   }
   bool isConsistent() const
   {
      if( m_rownames.num() != nRows() )
         return MSGinconsistent("SPxSIP");

      if( m_colnames.num() != nCols() )
         return MSGinconsistent("SPxSIP");      

      return SoPlex::isConsistent();
   }
};

SPxSIP* SPxSIP::fromSIPLP(SIPLP lptr)
{
   SPxSIP* spx = reinterpret_cast<SPxSIP*>(lptr); /*lint !e740*/

   assert(spx->isConsistent());

   return spx;
}

extern "C" int SIPstrongbranch(
   FILE*     ferr, 
   SIPInfaLP infaLP, 
   SIPLP     lptr, 
   double*   psol,
   int*      cand, 
   int       ncand,
   double*   down, 
   double*   up, 
   int       itlim)
{
   METHOD("SIPstrongbranch");

   assert(ferr   != NULL);
   assert(infaLP != NULL);
   assert(lptr   != NULL);
   assert(psol   != NULL);
   assert(cand   != NULL);
   assert(ncand  >  0);
   assert(down   != NULL);
   assert(up     != NULL);
   assert(itlim  >  0);

   SPxSIP* spx = SPxSIP::fromSIPLP(lptr);

   SoPlex::VarStatus* rowstat  = new SoPlex::VarStatus[spx->nRows()];
   SoPlex::VarStatus* colstat  = new SoPlex::VarStatus[spx->nCols()];
   double             oldBound;
   int                oldItlim = spx->terminationIter();
   int                c;
   SoPlex::Status     status = SoPlex::UNKNOWN;
   bool               error = false;

   spx->getBasis( rowstat, colstat );
   spx->setTerminationIter( itlim );
   spx->setType( SoPlex::LEAVE );
   
   for( c = 0; c < ncand && !error; ++c )
   {
      assert( SIP_FRAC( psol[cand[c]] ) );

      /* down branch */
      oldBound = spx->upper( cand[c] );
      spx->changeUpper( cand[c], SIP_DOWN( psol[cand[c]] ) );

#ifdef DEBUG_FILEOUTPUT
      spx->dumpFile( "strongdownLP.lp" );
#endif

      status = spx->solve();
      switch( status )
      {
      case SoPlex::ABORT_TIME:
      case SoPlex::ABORT_ITER:
      case SoPlex::OPTIMAL:
         down[c] = spx->value();
         break;
      case SoPlex::ABORT_VALUE:
      case SoPlex::INFEASIBLE:
         down[c] = spx->terminationValue();
         break;
      default:
         error = true;
         break;
      }
      spx->changeUpper( cand[c], oldBound );
      spx->setBasis( rowstat, colstat );

      if( !error )
      {
         /* up branch */
         oldBound = spx->lower( cand[c] );
         spx->changeLower( cand[c], SIP_CEIL( psol[cand[c]] ) );

#ifdef DEBUG_FILEOUTPUT
      spx->dumpFile( "strongupLP.lp" );
#endif

         status = spx->solve();
         switch( status )
         {
         case SoPlex::ABORT_TIME:
         case SoPlex::ABORT_ITER:
         case SoPlex::OPTIMAL:
            up[c] = spx->value();
            break;
         case SoPlex::ABORT_VALUE:
         case SoPlex::INFEASIBLE:
            up[c] = spx->terminationValue();
            break;
         case SoPlex::UNBOUNDED:
         default:
            error = true;
            break;
         }
         spx->changeLower( cand[c], oldBound );
         spx->setBasis( rowstat, colstat );
      }
   }

   spx->setTerminationIter( oldItlim );

   delete [] rowstat;
   delete [] colstat;

   if( error )
   {
      fprintf( ferr, "Error in SIPstrongbranch():" );
      fprintf( ferr, " SOPLEX status = %d\n", int(status) );
      return SIP_LPERROR;
   }
   return SIP_OKAY;
}

extern "C" int SIPoptLP(
   SIPInfaLP infaLP, 
   SIPLP     lptr)
{
   METHOD("SIPoptLP");

   assert(infaLP != NULL);
   assert(lptr   != NULL);

   SPxSIP* spx = SPxSIP::fromSIPLP(lptr);

   /* spx->setType(SoPlex::LEAVE); */
   
#if defined(DEBUG_FILEOUTPUT)
   spx->dumpFile( "subproblemLP.lp" );
#endif

   SoPlex::Status status = spx->solve();

   switch( status )
      {
      case SoPlex::ABORT_TIME:
      case SoPlex::ABORT_ITER:
      case SoPlex::ABORT_VALUE:
      case SoPlex::SINGULAR:
      case SoPlex::REGULAR:
      case SoPlex::UNKNOWN:
      case SoPlex::OPTIMAL:
      case SoPlex::UNBOUNDED:
      case SoPlex::INFEASIBLE:
         return SIP_OKAY;
      case SoPlex::NO_PROBLEM:
      case SoPlex::RUNNING:
      case SoPlex::ERROR:
      default:
         return SIP_LPERROR;
      }
}

extern "C" int SIPoptLPPrimal(
   SIPInfaLP infaLP, 
   SIPLP     lptr)
{
   METHOD("SIPoptLPPrimal");

   assert(infaLP != NULL);
   assert(lptr   != NULL);

   SPxSIP* spx = SPxSIP::fromSIPLP(lptr);

   /* spx->setType(SoPlex::ENTER); */
   
#ifdef DEBUG_FILEOUTPUT
   spx->dumpFile( "subproblemLP.lp" );
#endif

   SoPlex::Status status = spx->solve();

   switch( status )
      {
      case SoPlex::ABORT_TIME:
      case SoPlex::ABORT_ITER:
      case SoPlex::ABORT_VALUE:
      case SoPlex::SINGULAR:
      case SoPlex::REGULAR:
      case SoPlex::UNKNOWN:
      case SoPlex::OPTIMAL:
      case SoPlex::UNBOUNDED:
      case SoPlex::INFEASIBLE:
         return SIP_OKAY;
      case SoPlex::NO_PROBLEM:
      case SoPlex::RUNNING:
      case SoPlex::ERROR:
      default:
         return SIP_LPERROR;
      }
}

extern "C" int SIPgetLPStatus(
   SIPInfaLP infaLP,
   SIPLP     lptr,
   int       *solstat)
{
   METHOD("SIPgetLPStatus");

   assert(infaLP != NULL);
   assert(lptr   != NULL);

   SPxSIP* spx = SPxSIP::fromSIPLP(lptr);

   *solstat = int(spx->getStatus());
   return SIP_OKAY;
}


/** returns TRUE iff error occured during LP solve */
extern "C" int SIPerrorLP(
   SIPLP lptr,
   int solstat)
{
   METHOD("SIPerrorLP");

   assert(lptr != NULL);

   return (solstat == int(SoPlex::ERROR))
      ||  (solstat == int(SoPlex::NO_PROBLEM))
      ||  (solstat == int(SoPlex::RUNNING))
      ||  (solstat == int(SoPlex::UNKNOWN));
}

/** returns TRUE iff LP solution is stable, ie no numerical troubles occured */
extern "C" int SIPisStable(
   SIPLP lptr,
   int solstat)
{
   METHOD("SIPisStable");

   assert(lptr != NULL);

   return (solstat != int(SoPlex::SINGULAR));
}

/** returns TRUE iff LP is primal infeasible */
extern "C" int SIPisPrimalInfeasible(
   SIPLP lptr,
   int solstat)
{
   METHOD("SIPisPrimalInfeasible");

   assert(lptr != NULL);

   return (solstat == int(SoPlex::INFEASIBLE));
}

/** returns TRUE iff LP is primal unbounded */
extern "C" int SIPisPrimalUnbounded(
   SIPLP lptr,
   int solstat)
{
   METHOD("SIPisPrimalUnbounded");

   assert(lptr != NULL);

   return (solstat == int(SoPlex::UNBOUNDED));
}

/** returns TRUE iff LP is solved to optimality */
extern "C" int SIPisOptimal(
   SIPLP lptr,
   int solstat)
{
   METHOD("SIPisOptimal");

   assert(lptr != NULL);

   return (solstat == int(SoPlex::OPTIMAL));
}

/** returns TRUE iff a dual feasible solution has been found */
extern "C" int SIPisDualValid(
   SIPLP lptr,
   int solstat)
{
   METHOD("SIPisDualValid");

   assert(lptr != NULL);

   return (solstat == int(SoPlex::OPTIMAL));
}

/** returns TRUE iff objective limit is exceeded */
extern "C" int SIPexObjlim(
   SIPLP lptr,
   int solstat)
{
   METHOD("SIPexObjlim");

   assert(lptr != NULL);

   return (solstat == int(SoPlex::UNBOUNDED))
      ||  (solstat == int(SoPlex::ABORT_VALUE));
}

/** TRUE iff iteration limit has been exceeded */
extern "C" int SIPiterlim(
   SIPLP lptr,
   int solstat)
{
   METHOD("SIPiterlim");

   assert(lptr != NULL);

   return (solstat == int(SoPlex::ABORT_ITER));
}

/** TRUE iff time limit has been exceeded */
extern "C" int SIPtimelim(
   SIPLP lptr,
   int solstat)
{
   METHOD("SIPtimelim");

   assert(lptr != NULL);

   return (solstat == int(SoPlex::ABORT_TIME));
}


extern "C" int SIPopenInfa(
   FILE*      /*ferr*/, 
   SIPInfaLP* infaLP)
{
   METHOD("SIPopenInfa");

   assert(infaLP != NULL);

   *infaLP = reinterpret_cast<SIPInfaLP>(&dummy_environment); /*lint !e740*/
   
   return SIP_OKAY;
      
}

extern "C" int SIPfreeInfa(
   FILE*      /*ferr*/, 
   SIPInfaLP* infaLP)
{
   METHOD("SIPfreeInfa");

   assert(infaLP  != NULL);
   assert(*infaLP != NULL);

   *infaLP = NULL;

   return SIP_OKAY;
   
}

extern "C" int SIPopenLP(
   FILE*     /*ferr*/, 
   SIPInfaLP infaLP, 
   SIPLP*    lptr, 
   char*     name)
{
   METHOD("SIPopenLP");

   assert(infaLP != NULL);
   assert(lptr   != NULL);

   SPxSIP* spx = new SPxSIP;
   
   if (name != NULL)
      spx->setProbname(name);

   *lptr = reinterpret_cast<SIPLP>(spx);  /*lint !e740*/

   return SIP_OKAY;

}

extern "C" int SIPfreeLP(
   FILE*     /*ferr*/, 
   SIPInfaLP infaLP, 
   SIPLP*    lptr)
{
   METHOD("SIPfreeLP");

   assert(infaLP != NULL);
   assert(lptr   != NULL);
   
   if (*lptr != NULL)
   {
      delete SPxSIP::fromSIPLP(*lptr);  /*lint !e740*/
      *lptr = NULL;
   }

   return SIP_OKAY;
   
}

extern "C" int SIPcopyLP(
   FILE*     ferr, 
   SIPInfaLP infaLP, 
   SIPLP     lptr, 
   int       ncol, 
   int       nrow,
   int       objsen, 
   double*   obj, 
   double*   rhs, 
   char*     sen, 
   int*      beg,
   int*      cnt, 
   int*      ind, 
   double*   val, 
   double*   lb, 
   double*   ub,
   char**    cname, 
   char**    rname)
{
   char      defname[128];

   METHOD("SIPcopyLP");

   assert(ferr   != NULL);
   assert(infaLP != NULL);
   assert(lptr   != NULL);
   assert(ncol   >  0);
   assert(nrow   >  0);
   assert(objsen == -1 || objsen == +1);
   assert(obj    != NULL);
   assert(rhs    != NULL);
   assert(sen    != NULL);
   assert(beg    != NULL);
   assert(cnt    != NULL);
   assert(ind    != NULL);
   assert(val    != NULL);
   assert(lb     != NULL);
   assert(ub     != NULL);
   /*   assert(cname  != NULL);*/
   /*   assert(rname  != NULL);*/
 
   SPxSIP*   spx = SPxSIP::fromSIPLP(lptr);
   LPColSet  cols(ncol);
   LPRowSet  rows(nrow);
   NameSet   colnames(ncol);
   NameSet   rownames(nrow);
   DSVector  colVector(nrow);
   DSVector  emptyVector(0);
   int       i;

   spx->clearWithNames();

   for(i = 0; i < nrow; i++)
   {
      switch(sen[i])
      {
      case 'L':
         rows.add(-soplex::infinity, emptyVector, rhs[i]);
         break;
      case 'G':
         rows.add(rhs[i], emptyVector, soplex::infinity);
         break;
      case 'E':
         rows.add(rhs[i], emptyVector, rhs[i]);
         break;
      case 'R':
#if 0
         assert(rng != NULL);
         if (rng[i] > 0)
            rows.add(rhs[i], emptyVector, rhs[i] + rng[i]);
         else
            rows.add(rhs[i] + rng[i], emptyVector, rhs[i]);
         break;
#else
         std::cout << "Ranges are not supported!" << std::endl;
         abort();
#endif
      default:
         abort();
      }
      if( rname != NULL )
         rownames.add(rname[i]);
      else
      {
         sprintf( defname, "r%d", i );
         rownames.add( defname );
      }
   }
   spx->addRowsWithNames(rows, rownames);

   for(i = 0; i < ncol; i++)
   {
      colVector.clear();

      for(int j = 0; j < cnt[i]; j++)
         colVector.add(ind[beg[i] + j], val[beg[i] + j]);

      cols.add(obj[i], lb[i], colVector, ub[i]);
      if( cname != NULL )
         colnames.add(cname[i]);
      else
      {
         sprintf( defname, "x%d", i );
         colnames.add( defname );
      }
   }
   spx->changeSense(objsen == 1 ? SPxLP::MINIMIZE : SPxLP::MAXIMIZE);
   spx->addColsWithNames(cols, colnames);
   
   /* we do not use the names for now. They where only needed for the
    * output routines, and those should be part of SIP anyway.
    */

#if DEBUGGING >= 1
   std::cerr << *spx;
#endif
#if defined(DEBUG_FILEOUTPUT)
   spx->dumpFile( "copiedLP.lp" );
#endif

   return SIP_OKAY;
   
}

extern "C" int SIPsetbase(
   FILE*     /*ferr*/, 
   SIPInfaLP infaLP, 
   SIPLP     lptr, 
   double*   /*dnorm*/,
   int*      cstat, 
   int*      rstat, 
   int       /*pricing*/)
{
   METHOD("SIPsetbase");

   int i;

   assert(infaLP != NULL);
   assert(lptr   != NULL);
   assert(cstat  != NULL);
   assert(rstat  != NULL);

   SPxSIP* spx = SPxSIP::fromSIPLP(lptr); 
   SoPlex::VarStatus* spxcstat = new SoPlex::VarStatus[spx->nCols()];
   SoPlex::VarStatus* spxrstat = new SoPlex::VarStatus[spx->nRows()];

#if DEBUGGING >= 2
   print_base('<', spx->nCols(), spx->nRows(), cstat, rstat);
#endif

   for(i = 0; i < spx->nRows(); i++)
   {
      switch(rstat[i])
      {
      case BASE_AT_LOWER:
         spxrstat[i] = SoPlex::ON_LOWER;
         break;
      case BASE_IS_BASIC:
         spxrstat[i] = SoPlex::BASIC;
         break;
      case BASE_AT_UPPER:
         spxrstat[i] = SoPlex::ON_UPPER;
         break;
      case BASE_IS_FREE:
         spxrstat[i] = SoPlex::ZERO;
         break;
      default:
         abort();
      }
   }

   for(i = 0; i < spx->nCols(); i++)
   {
      switch(cstat[i])
      {
      case BASE_AT_LOWER:
         spxcstat[i] = SoPlex::ON_LOWER;
         break;
      case BASE_IS_BASIC:
         spxcstat[i] = SoPlex::BASIC;
         break;
      case BASE_AT_UPPER:
         spxcstat[i] = SoPlex::ON_UPPER;
         break;
      case BASE_IS_FREE:
         spxcstat[i] = SoPlex::ZERO;
         break;
      default:
         abort();
      }
   }
   spx->setBasis(spxrstat, spxcstat);

   delete[] spxcstat;
   delete[] spxrstat;
   
   return SIP_OKAY;
}

extern "C" int SIPgetbase(
   FILE*     /*ferr*/, 
   SIPInfaLP infaLP, 
   SIPLP     lptr,
   double*   /*dnorm*/, 
   int*      cstat, 
   int*      rstat)
{
   METHOD("SIPgetbase");

   assert(infaLP != NULL);
   assert(lptr   != NULL);

   SPxSIP* spx = SPxSIP::fromSIPLP(lptr);

   if (rstat != NULL)
   {
      for(int i = 0; i < spx->nRows(); i++)
      {
         switch(spx->getBasisRowStatus(i))
         {
         case SoPlex::BASIC:
            rstat[i] = BASE_IS_BASIC;
            break;	  
         case SoPlex::FIXED:
         case SoPlex::ON_LOWER:
            rstat[i] = BASE_AT_LOWER;
            break;
         case SoPlex::ON_UPPER:
            rstat[i] = BASE_AT_UPPER;
            break;
         case SoPlex::ZERO:
            rstat[i] = BASE_IS_FREE;
            break;
         default:
            abort();
         }
      }
   }

   if (cstat != NULL)
   {
      for(int i = 0; i < spx->nCols(); i++)
      {
         switch(spx->getBasisColStatus(i))
         {
         case SoPlex::BASIC:
            cstat[i] = BASE_IS_BASIC;
            break;	  
         case SoPlex::FIXED:
         case SoPlex::ON_LOWER:
            cstat[i] = BASE_AT_LOWER;
            break;
         case SoPlex::ON_UPPER:
            cstat[i] = BASE_AT_UPPER;
            break;
         case SoPlex::ZERO:
            cstat[i] = BASE_IS_FREE;
            break;
         default:
            abort();
         }
      }
   }
#if DEBUGGING >= 2
   print_base('>', spx->nCols(), spx->nRows(), cstat, rstat);
#endif

   return SIP_OKAY;
}

extern "C" int SIPgetsol(
   FILE*     /*ferr*/, 
   SIPInfaLP infaLP, 
   SIPLP     lptr, 
   double*   objval,
   double*   psol, 
   double*   pi, 
   double*   slck, 
   double*   redcost)
{
   METHOD("SIPgetsol");

   assert(infaLP != NULL);
   assert(lptr   != NULL);
    
   SPxSIP* spx = SPxSIP::fromSIPLP(lptr);

   if (objval != NULL)
      *objval = spx->value();

   if (psol != NULL)
   {
      Vector tmp(spx->nCols(), psol);
      spx->getPrimal(tmp);
   }
   if (redcost != NULL)
   {
      Vector tmp(spx->nCols(), redcost);
      spx->getRdCost(tmp);
   }
   if ((slck != NULL) || (pi != NULL))
   {
      int     rows  = spx->nRows() ;
      DVector tmp(rows);

      if (slck != NULL)
      {
         spx->getSlacks(tmp);  /* returns  tmp = Ax */

         for(int i = 0; i < rows; i++)
         {
            double lhs = spx->lhs(i);  /* SOPLEX: lhs <= act(x) <= rhs */
            double rhs = spx->rhs(i);
            double act = tmp[i];
            
            if( lhs <= -soplex::infinity || /*  <= row  -> slack(x) = rhs - act(x)  in  [0,inf] */
               lhs == rhs               )   /*  == row  -> slack(x) = rhs - act(x)  in  [0,0] */
               slck[i] = rhs - act;
            else                            /*  >= row  -> slack(x) = act(x) - lhs  in  [0,inf] */
               slck[i] = act - lhs;         /* RNG row  -> slack(x) = act(x) - lhs  in  [0,rng] */
         }
      }
      if (pi != NULL)
      {
         spx->getDual(tmp);
         for(int i = 1; i < rows; i++)
            pi[i - 1] = tmp[i];
      }
   }
   return SIP_OKAY;

}

extern "C" void SIPsetintparLP(
   SIPInfaLP infaLP, 
   SIPLP     lptr,
   int       type, 
   int       ival)
{
   METHOD("SIPsetintparLP");
#if DEBUGGING >= 3
   std::cerr << "Type " << type << " " << ival << std::endl;
#endif

   assert(infaLP != NULL);
   assert(lptr != NULL);

   SPxSIP* spx = SPxSIP::fromSIPLP(lptr);

   switch(type) 
   {
   case SIP_FROMSCRATCH:
      spx->setFromScratch(ival == SIP_ON);
      break;
   case SIP_LPITLIM:
      spx->setTerminationIter(ival);
      break;
   case SIP_FASTMIP:
      break;
   case SIP_PRICING:
      break;
   case SIP_LPINFO:
      break;
   default:
      abort();
   }
}

extern "C" void SIPsetdblparLP(
   SIPInfaLP infaLP, 
   SIPLP     lptr,
   int       type, 
   double    dval)
{
   METHOD("SIPsetdblparLP");
#if DEBUGGING >= 3
   std::cerr << "Type " << type << " " << dval << std::endl;
#endif

   assert(infaLP != NULL);
   assert(lptr != NULL);

   SPxSIP* spx = SPxSIP::fromSIPLP(lptr);

   switch(type) 
   {
   case SIP_FEASTOL:
      spx->setDelta(dval);
      break;
   case SIP_LOBJLIM:
      spx->setObjLoLimit(dval);
      break;
   case SIP_UOBJLIM:
      spx->setObjUpLimit(dval);
      /* if (dval == SIP_INFINITY)
       *    spx->setTerminationValue(SoPlex::infinity);
       * else
       *    spx->setTerminationValue(dval);
       */
      break;
   case SIP_LPTILIM:
      spx->setTerminationTime(dval);
      break;
   default:
      abort();
   }
}

extern "C" int SIPgetintparLP(
   SIPInfaLP infaLP, 
   SIPLP     lptr,
   int       type, 
   int*      ival)
{
   METHOD("SIPgetintparLP");

   assert(infaLP != NULL);
   assert(lptr != NULL);
   assert(ival != NULL);

   SPxSIP* spx = SPxSIP::fromSIPLP(lptr);
   int restat  = SIP_OKAY;

   switch(type) 
   {
   case SIP_FROMSCRATCH:
      *ival = spx->getFromScratch() ? SIP_ON : SIP_OFF;
      break;
   case SIP_LPNROW:      
      *ival = spx->nRows();
      break;
   case SIP_LPIT:
      *ival = spx->basis().iteration();
      break;
   case SIP_LPIT1:
      /* There is no phase I in SoPlex */
      *ival = 0;
      break;
   default:
      restat = SIP_LPERROR;
      break;
   }
   return restat;
}

extern "C" int SIPgetdblparLP(
   SIPInfaLP infaLP, 
   SIPLP     lptr,
   int       type, 
   double*   dval)
{
   METHOD("SIPgetdblparLP");

   assert(infaLP != NULL);
   assert(lptr != NULL);
   assert(dval != NULL);

   SPxSIP* spx = SPxSIP::fromSIPLP(lptr);
   int restat  = SIP_OKAY;

   switch(type)
   {
   case SIP_FEASTOL:
      *dval = spx->delta();
      break;
   default:
      restat = SIP_LPERROR;
      break;
   }
   return restat;
}

extern "C" int SIPgetDefaultLPFeastol(
   SIPInfaLP infaLP,
   double* LPFeastol)
{
   assert(infaLP != NULL);
   *LPFeastol = DEFAULT_BND_VIOL;
   return SIP_OKAY;
}


extern "C" int SIPwriteLP(
   FILE*     /*ferr*/, 
   SIPInfaLP infaLP, 
   SIPLP     lptr, 
   char*     fname)
{
   METHOD("SIPwriteLP");

   assert(infaLP != NULL);
   assert(lptr   != NULL);
   assert(fname  != NULL);

   SPxSIP* spx = SPxSIP::fromSIPLP(lptr);

   spx->dumpFile(fname);

   return SIP_OKAY;
}

extern "C" int SIPwriteB(
   FILE*     /*ferr*/, 
   SIPInfaLP infaLP, 
   SIPLP     lptr, 
   char*     fname)
{
   METHOD("SIPwriteB");

   assert(infaLP != NULL);
   assert(lptr   != NULL);
   assert(fname  != NULL);

   SPxSIP* spx = SPxSIP::fromSIPLP(lptr);

   spx->dumpBasis(fname);

   return SIP_OKAY;
}


extern "C" int SIPgetlb(
   FILE*     /*ferr*/, 
   SIPInfaLP infaLP, 
   SIPLP     lptr,
   double*   lb, 
   int       beg, 
   int       end)
{
   METHOD("SIPgetlb");

   assert(infaLP != NULL);
   assert(lptr   != NULL);
   assert(lb     != NULL);

   SPxSIP* spx = SPxSIP::fromSIPLP(lptr);

   assert(beg    >= 0);
   assert(end    >= beg);
   assert(end    <  spx->nCols());

   for(int i = beg; i <= end; i++)
      lb[i - beg] = spx->lower(i);

   return SIP_OKAY; 
}

extern "C" int SIPgetub(
   FILE*     /*ferr*/, 
   SIPInfaLP infaLP, 
   SIPLP     lptr,
   double*   ub, 
   int       beg, 
   int       end)
{
   METHOD("SIPgetub");

   assert(infaLP != NULL);
   assert(lptr   != NULL);
   assert(ub     != NULL);

   SPxSIP* spx = SPxSIP::fromSIPLP(lptr);

   assert(beg    >= 0);
   assert(end    >= beg);
   assert(end    <  spx->nCols());
 
   for(int i = beg; i <= end; i++)
      ub[i - beg] = spx->upper(i);

   return SIP_OKAY;
}

extern "C" int SIPchgbds(
   FILE*     ferr, 
   SIPInfaLP infaLP, 
   SIPLP     lptr,
   int       cnt, 
   int*      idx, 
   char*     lu, 
   double*   bd)
{
   METHOD("SIPchgbds");

   assert(infaLP != NULL);
   assert(lptr   != NULL);
   assert(cnt    >  0);
   assert(idx    != NULL);
   assert(lu     != NULL);
   assert(bd     != NULL);

   SPxSIP* spx = SPxSIP::fromSIPLP(lptr);

   for(int i = 0; i < cnt; i++)
   {
      assert(idx[i] >= 0);
      assert(idx[i] <  spx->nCols());

      switch( lu[i] )
      {
      case 'U':
         spx->changeUpper(idx[i], bd[i]);
         break ;
      case 'L':
         spx->changeLower(idx[i], bd[i]);
         break ;
      case 'B':
         spx->changeUpper(idx[i], bd[i]);
         spx->changeLower(idx[i], bd[i]);
         break ;
      default:
         fprintf( ferr, "Error in SIPchgbds(): " );
         fprintf( ferr, "Unknown bound symbol <%c>.\n", lu[i] );
         abort();
      }
   }
   return SIP_OKAY;
}

extern "C" int SIPchgrhs(
   FILE*     /*ferr*/, 
   SIPInfaLP infaLP, 
   SIPLP     lptr,
   int       cnt, 
   int*      idx, 
   double*   rhs)
{
   METHOD("SIPchgbds");

   assert(infaLP != NULL);
   assert(lptr   != NULL);
   assert(cnt    >  0);
   assert(idx    != NULL);
   assert(rhs    != NULL);

   SPxSIP* spx = SPxSIP::fromSIPLP(lptr);

   for(int i = 0; i < cnt; i++)
   {
      assert(idx[i] >= 0);
      assert(idx[i] <  spx->nRows());
      spx->changeRhs(idx[i], rhs[i]);
   }
   return SIP_OKAY;
}

extern "C" void SIPchgobjsen(
   FILE*     /*ferr*/, 
   SIPInfaLP infaLP, 
   SIPLP     lptr, 
   int       objsen)
{
   METHOD("SIPchgobjsen");

   assert(infaLP != NULL);
   assert(lptr   != NULL);
   assert(objsen == -1 || objsen == +1);

   SPxSIP* spx = SPxSIP::fromSIPLP(lptr);

   spx->changeSense(objsen == 1 ? SoPlex::MINIMIZE : SoPlex::MAXIMIZE);
}

extern "C" int SIPdelrows(
   FILE*     /*ferr*/, 
   SIPInfaLP infaLP, 
   SIPLP     lptr, 
   int*      dstat)
{
   METHOD("SIPdelrows");

   assert(infaLP != NULL);
   assert(lptr   != NULL);
   assert(dstat  != NULL);

   SPxSIP* spx = SPxSIP::fromSIPLP(lptr);

   spx->removeRowsWithNames( dstat );

   return SIP_OKAY;
}


/* returns the indices of the basis variables.
 * head[b] < 0 => basis variable b is the slack variable
                  belonging to row -head[b]-1.
 */
extern "C" int SIPgetBind(
   FILE*     /*ferr*/, 
   SIPInfaLP infaLP, 
   SIPLP     lptr,
   int*      head )
{
   METHOD("SIPgetBind");

   assert(infaLP != NULL);
   assert(head != NULL);

   SPxSIP* spx = SPxSIP::fromSIPLP(lptr);

   int i;

   for( i = 0; i < spx->nRows(); ++i )
      {
         SPxId id = spx->basis().baseId( i );
         if( spx->isId( id ) ) /* column id? */
            head[i] = spx->number( id );
         else                  /* row id?    */
            head[i] = - spx->number( id ) - 1;
      }

   return SIP_OKAY;
}

/* returns a row of the constraint matrix A */
extern "C" int SIPgetrow(
   FILE*     ferr, 
   SIPInfaLP infaLP, 
   SIPLP     lptr,
   int       i, 
   double*   val, 
   int*      ind, 
   int*      nnonz)
{
   METHOD("SIPgetrow");

   assert(ferr   != NULL);
   assert(infaLP != NULL);
   assert(lptr   != NULL);
   assert(val    != NULL);
   assert(ind    != NULL);
   assert(nnonz  != NULL);

   SPxSIP* spx = SPxSIP::fromSIPLP(lptr);

   assert(i     >= 0);
   assert(i     <  spx->nRows());

   /* needed because of creation of Svector */
   {
      DSVector  row( spx->rowVector(i) );
      int       cnt = 0;     
      
      /* sort row vector, because the pool ieqs are sorted, and
         SIPcheckpoollp fails if SIPgetrow returns an unsorted
         row vector */
      row.sort();

      assert(*nnonz >= row.size());
      
      for(int j = 0; j < row.size(); j++)
         {
            ind[cnt] = row.index(j);
            val[cnt] = row.value(j);
            cnt++;
         }
      *nnonz = cnt;
   }
   return SIP_OKAY;
}

/** returns a row of the inverse B^-1 of the basis matrix B */
extern "C" int SIPgetrowBinv(
   FILE*     ferr, 
   SIPInfaLP infaLP, 
   SIPLP     lptr,
   int       i, 
   double*   val )
{
   METHOD("SIPgetrowBinv");

   assert(ferr   != NULL);
   assert(infaLP != NULL);
   assert(lptr   != NULL);
   assert(val    != NULL);

   SPxSIP* spx = SPxSIP::fromSIPLP(lptr);

   assert(i     >= 0);
   assert(i     <  spx->nRows());
      
   /* solve system "x = e_i^T * B^-1" to get i'th row of B^-1 */
   soplex::Vector x( spx->nRows(), val );
   spx->basis().coSolve( x, spx->unitVector( i ) );

   return SIP_OKAY;
}


/** returns a row of B^-1 * A */
extern "C" int SIPgetrowBinvA(
   FILE*     ferr, 
   SIPInfaLP infaLP, 
   SIPLP     lptr,
   int       i, 
   double*   binv,
   double*   val )
{
   int restat = SIP_OKAY;
   int freebinv = FALSE;
   int col;

   METHOD("SIPgetrowBinvA");

   assert(ferr   != NULL);
   assert(infaLP != NULL);
   assert(lptr   != NULL);
   assert(val    != NULL);

   SPxSIP* spx = SPxSIP::fromSIPLP(lptr);

   assert(i >= 0);
   assert(i <  spx->nRows());

   if( binv == NULL )
   {
      ALLOC_OKAY( binv = static_cast<double*>(allocMemoryArrayCPP(spx->nRows(), sizeof(*binv))) );
      freebinv = TRUE;
      restat = SIPgetrowBinv( ferr, infaLP, lptr, i, binv );
      if( restat != SIP_OKAY )
         goto TERMINATE;
   }

   assert( binv != NULL );      

   {
      soplex::Vector binvvec( spx->nRows(), binv );
      
      for( col = 0; col < spx->nCols(); ++col )
         val[col] = binvvec * spx->colVector( col );  /* scalar product */
   }

 TERMINATE:   
   if( freebinv )
      freeMemoryArray( binv );

   return restat;
}


extern "C" int SIPaddrow(
   FILE*     ferr, 
   SIPInfaLP infaLP, 
   SIPLP     lptr, 
   int       rcnt, 
   int       nzcnt,
   double*   rhs, 
   char*     sns, 
   int*      beg, 
   int*      ind, 
   double*   val,
   char**    name)
{
   METHOD("SIPaddrow");

   assert(ferr   != NULL);
   assert(infaLP != NULL);
   assert(lptr   != NULL);
   assert(rhs    != NULL);
   assert(sns    != NULL);
   assert(beg    != NULL);
   assert(ind    != NULL);
   assert(val    != NULL);
   assert(name   != NULL);
   assert(rcnt   >= 0);
   assert(nzcnt  >= 0);

   SPxSIP* spx = SPxSIP::fromSIPLP(lptr);

   DSVector rowvec(spx->nCols());
   LPRowSet rset(rcnt, nzcnt);
   NameSet  rownames(rcnt);

   for(int i = 0; i < rcnt; i++)
   {
      rowvec.clear();

      int end = (i < rcnt - 1) ? beg[i + 1] : nzcnt;
      
      for(int k = beg[i]; k < end; k++)
         rowvec.add(ind[k], val[k]);

      switch(sns[i])
      {
      case 'L':
         rset.add(-soplex::infinity, rowvec, rhs[i]);
         rownames.add(name[i]);
         break;
      case 'G':
         rset.add(rhs[i], rowvec, soplex::infinity);
         rownames.add(name[i]);
         break;
      case 'E':
         rset.add(rhs[i], rowvec, rhs[i]);
         rownames.add(name[i]);
         break;
      case 'R':
         /* ranged rows not supported */
      default:
         abort();
      }
   }   
   spx->addRowsWithNames(rset, rownames);

   return SIP_OKAY;
}

extern "C" int SIPreadFile(
   SET *set,
   MIP *mip,
   const char *filename)
{
   SPxLP         lp;
   NameSet       colnames;
   NameSet       rownames;
   DIdxSet       intvars;
   int           i;
   int           k;
   
   METHOD("SIPreadFile");

   assert( set != NULL );
   assert( mip != NULL );
   assert( filename != NULL );
   
   /* check filename */
   if( strlen(filename) == 0 )
   {
      fprintf (stderr, "No file specified\n");
      return SIP_NOFILE;
   }

   /* store filename in MIP structure */
   ALLOC_OKAY( mip->name = static_cast<char*>(allocMemoryArrayCPP(strlen(filename)+1, sizeof(*mip->name))) );
   strcpy( mip->name, filename );

   /* read file into SPxLP object */
   std::ifstream file(filename);

   if (!file || !lp.read(file, &rownames, &colnames, &intvars))
   {
      fprintf( stderr, "Could not read file %s\n", filename);
      return SIP_NOFILE;
   }
   rownames.memPack();
   colnames.memPack();

#if DEBUGGING >= 1
   std::cerr << "As read by SoPlex::spxlpfread" << std::endl;
   std::cerr << lp;
#endif

   /* get matrix size */
   int numnz = 0;
   for(i = 0; i < lp.nCols(); i++)
      numnz += lp.colVector(i).size();

   /* get number of columns */
   mip->ncol   = lp.nCols();

   /* get matrix in column form */
   ALLOC_OKAY( mip->col = static_cast<MAT*>(allocMemoryCPP(sizeof(*mip->col))) );
   ALLOC_OKAY( mip->col->pnt = static_cast<int*>(allocMemoryArrayCPP(mip->ncol+1, sizeof(*mip->col->pnt))) );
   ALLOC_OKAY( mip->col->ind = static_cast<int*>(allocMemoryArrayCPP(numnz, sizeof(*mip->col->ind))) );
   ALLOC_OKAY( mip->col->val = static_cast<double*>(allocMemoryArrayCPP(numnz, sizeof(*mip->col->val))) );

   for(k = 0, i = 0; i < mip->ncol; i++)
   {
      SVector col = lp.colVector(i);
      
      mip->col->pnt[i] = k;

      for(int j = 0; j < col.size(); j++)
      {
         mip->col->ind[k] = col.index(j);
         mip->col->val[k] = col.value(j);
         k++;
      }
   }
   mip->col->pnt[mip->ncol] = numnz;

   assert(k == numnz);

   /* get objective function */
   ALLOC_OKAY( mip->obj = static_cast<double*>(allocMemoryArrayCPP(mip->ncol, sizeof(*mip->obj))) );

   /* assign memory position to tmpobj vector */
   Vector tmpobj(mip->ncol, mip->obj);  
   /* get objective function */
   lp.getObj( tmpobj );                 

   /* get variable bounds */
   ALLOC_OKAY( mip->lb = static_cast<double*>(allocMemoryArrayCPP(mip->ncol, sizeof(*mip->lb))) );
   ALLOC_OKAY( mip->ub = static_cast<double*>(allocMemoryArrayCPP(mip->ncol, sizeof(*mip->ub))) );

   for(i = 0; i < mip->ncol; i++)
   {
      mip->lb[i] = lp.lower(i);
      mip->ub[i] = lp.upper(i);
   }

   /* get column type */
   ALLOC_OKAY( mip->ctype = static_cast<char*>(allocMemoryArrayCPP(mip->ncol, sizeof(*mip->ctype))) );

   for(i = 0; i < mip->ncol; i++)
      mip->ctype[i] = 'C';

   for(i = 0; i < intvars.size(); i++)
   {
      int idx = intvars.index(i);

      mip->ctype[idx] = 
         (lp.lower(idx) == 0.0 && lp.upper(idx) == 1.0) ? 'B' : 'I';
   }

   /* get column names */
   assert(colnames.num() == mip->ncol);

   mip->cstoresz = colnames.memSize();

   ALLOC_OKAY( mip->cstore = static_cast<char*>(allocMemoryArrayCPP(mip->cstoresz, sizeof(*mip->cstore))) );
   ALLOC_OKAY( mip->cname = static_cast<char**>(allocMemoryArrayCPP(mip->cstoresz, sizeof(*mip->cname))) );

   k = 0;

   for(i = 0; i < mip->ncol; i++)
   {
      mip->cname[i] = &(mip->cstore[k]);
      strcpy(mip->cname[i], colnames[i]);
      k += (int)strlen(colnames[i]) + 1;
      assert(k <= mip->cstoresz);
   }

   /* get number of rows */
   mip->nrow = lp.nRows();

   /* get matrix in row form */
   mip->row = NULL;  /* The unpresolved matrix is not needed in row form */

   /* get objective sense */
   mip->objsen = lp.spxSense() == SPxLP::MAXIMIZE ? -1 : 1;  

   /* get row rhs, range, and sense */
   ALLOC_OKAY( mip->rhs = static_cast<double*>(allocMemoryArrayCPP(mip->nrow, sizeof(*mip->rhs))) );
   ALLOC_OKAY( mip->rngval = static_cast<double*>(allocMemoryArrayCPP(mip->nrow, sizeof(*mip->rngval))) );
   ALLOC_OKAY( mip->sen = static_cast<char*>(allocMemoryArrayCPP(mip->nrow, sizeof(*mip->sen))) );

   for( i = 0; i < mip->nrow; i++ )
   {
      switch(lp.rowType(i))
      {
      case LPRow::LESS_EQUAL :
         mip->sen   [i] = 'L';
         mip->rhs   [i] = lp.rhs(i);
         mip->rngval[i] = 0.0;
         break; 
      case LPRow::EQUAL :
         mip->sen   [i] = 'E';
         mip->rhs   [i] = lp.rhs(i);
         mip->rngval[i] = 0.0;
         break; 
      case LPRow::GREATER_EQUAL :
         mip->sen   [i] = 'G';
         mip->rhs   [i] = lp.lhs(i);
         mip->rngval[i] = 0.0;
         break; 
      case LPRow::RANGE :
         mip->sen   [i] = 'R';
         mip->rhs   [i] = lp.lhs(i);
         mip->rngval[i] = lp.rhs(i) - lp.lhs(i);
         break; 
      default :
         abort();
      }
   }   
   /* get row names */
   assert(rownames.num() == mip->nrow);

   mip->rstoresz = rownames.memSize();
   mip->rstoresz = rownames.memSize();
   ALLOC_OKAY( mip->rstore = static_cast<char*>(allocMemoryArrayCPP(mip->rstoresz, sizeof(*mip->rstore))) );
   ALLOC_OKAY( mip->rname = static_cast<char**>(allocMemoryArrayCPP(mip->rstoresz, sizeof(*mip->rname))) );

   k = 0;

   for(i = 0; i < mip->nrow; i++)
   {
      mip->rname[i] = &(mip->rstore[k]);
      strcpy(mip->rname[i], rownames[i]);
      k += (int)strlen(rownames[i]) + 1;
      assert(k <= mip->rstoresz);
   }

#if defined(DEBUG_FILEOUTPUT)
   SIPwriteIP( set, mip, "originalIP.lp", true );
#endif

   return SIP_OKAY;
}


#endif
