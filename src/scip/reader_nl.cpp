/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2021 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not visit scipopt.org.         */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   reader_nl.cpp
 * @brief  AMPL .nl file reader
 * @author Stefan Vigerske
 *
 * The code (including some comments) for the AMPLProblemHandler class based on ModelingSystemAMPL.cpp,
 * the nl-reader of SHOT (https://github.com/coin-or/SHOT).
 *
 * The code for SOS reading is based on the AMPL/Bonmin interface (https://github.com/coin-or/Bonmin).
 *
 * For documentation on ampl::mp, see https://ampl.github.io.
 */

/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <string>
#include <vector>

#include "scip/reader_nl.h"
#include "scip/cons_linear.h"
#include "scip/cons_nonlinear.h"
#include "scip/cons_sos1.h"
#include "scip/cons_sos2.h"
#include "scip/expr_var.h"
#include "scip/expr_value.h"
#include "scip/expr_sum.h"
#include "scip/expr_product.h"
#include "scip/expr_pow.h"
#include "scip/expr_log.h"
#include "scip/expr_exp.h"
#include "scip/expr_sin.h"
#include "scip/expr_cos.h"
#include "scip/expr_abs.h"

// disable -Wshadow warnings for upcoming includes of AMPL/MP
// disable -Wimplicit-fallthrough as I don't want to maintain extra comments in AMPL/MP code to suppress these
#ifdef __GNUC__
#pragma GCC diagnostic ignored "-Wshadow"
#pragma GCC diagnostic ignored "-Wimplicit-fallthrough"
#endif

#include "mp/nl-reader.h"
#include "mp/sol.h"

#define READER_NAME             "nlreader"
#define READER_DESC             "AMPL .nl file reader"
#define READER_EXTENSION        "nl"

#define SCIP_CALL_THROW(x) \
   do                                                                                                   \
   {                                                                                                    \
      SCIP_RETCODE throw_retcode;                                                                       \
      if( ((throw_retcode) = (x)) != SCIP_OKAY )                                                        \
         throw std::logic_error("Error <" + std::to_string((int)throw_retcode) + "> in function call"); \
   }                                                                                                    \
   while( false )

/*
 * Data structures
 */

/** problem data */
struct SCIP_ProbData
{
   SCIP_VAR**            vars;               /**< variables */
   int                   nvars;              /**< number of variables */
};

/*
 * Local methods
 */

// forward declaration
static SCIP_DECL_PROBDELORIG(probdataDelOrigNl);

class AMPLProblemHandler : public mp::NullNLHandler<SCIP_EXPR*>
{
private:
   SCIP* scip;
   SCIP_PROBDATA* probdata;

   std::vector<SCIP_Real>  nlconslhss;
   std::vector<SCIP_Real>  nlconsrhss;
   std::vector<SCIP_EXPR*> nlconsexprs;
   std::vector<SCIP_CONS*> linconss;
   SCIP_EXPR* objexpr;

   std::vector<SCIP_EXPR*> varexprs;
   std::vector<SCIP_EXPR*> exprstorelease;

public:
   AMPLProblemHandler(
      SCIP*              scip_,
      const char*        filename
      )
   : scip(scip_),
     probdata(NULL),
     objexpr(NULL)
   {
      assert(scip != NULL);
      assert(filename != NULL);

      const char* filebasename = strrchr(const_cast<char*>(filename), '/');
      if( filebasename == NULL )
         filebasename = filename;
      else
         ++filebasename;

      SCIP_CALL_THROW( SCIPallocClearMemory(scip, &probdata) );

      /* initialize empty SCIP problem */
      SCIP_CALL_THROW( SCIPcreateProb(scip, filebasename, probdataDelOrigNl, NULL, NULL, NULL, NULL, NULL, probdata) );
   }

   ~AMPLProblemHandler()
   {
      // finalize() must have been called
      assert(varexprs.empty());
      assert(exprstorelease.empty());
   }

   void OnHeader(
      const mp::NLHeader& h
      )
   {
      char name[SCIP_MAXSTRLEN];
      SCIP_CONS* cons;

      assert(probdata->vars == NULL);

      probdata->nvars = h.num_vars;
      SCIP_CALL_THROW( SCIPallocBlockMemoryArray(scip, &probdata->vars, probdata->nvars) );

      varexprs.resize(h.num_nl_vars_in_both + h.num_nl_vars_in_cons + h.num_nl_vars_in_objs);

      // create variables
      // create variable expressions for nonlinear variables
      for( int i = 0; i < h.num_vars; ++i )
      {
         SCIP_VARTYPE vartype;
         // Nonlinear variables in both constraints and objective
         if( i < h.num_nl_vars_in_both - h.num_nl_integer_vars_in_both )
            vartype = SCIP_VARTYPE_CONTINUOUS;
         else if( i < h.num_nl_vars_in_both )
            vartype = SCIP_VARTYPE_INTEGER;
         // Nonlinear variables in constraints
         else if( i - h.num_nl_vars_in_both < h.num_nl_vars_in_cons - h.num_nl_integer_vars_in_cons )
            vartype = SCIP_VARTYPE_CONTINUOUS;
         else if( i - h.num_nl_vars_in_both < h.num_nl_vars_in_cons )
            vartype = SCIP_VARTYPE_INTEGER;
         // Nonlinear variables in objective
         else if( i - h.num_nl_vars_in_both - h.num_nl_vars_in_cons < h.num_nl_vars_in_objs - h.num_nl_integer_vars_in_objs )
            vartype = SCIP_VARTYPE_CONTINUOUS;
         else if( i - h.num_nl_vars_in_both - h.num_nl_vars_in_cons < h.num_nl_vars_in_objs )
            vartype = SCIP_VARTYPE_INTEGER;
         // Linear variables
         else if( i < h.num_vars - h.num_linear_binary_vars - h.num_linear_integer_vars )
            vartype = SCIP_VARTYPE_CONTINUOUS;
         else if( i < h.num_vars - h.num_linear_integer_vars )
            vartype = SCIP_VARTYPE_BINARY;
         else
            vartype = SCIP_VARTYPE_INTEGER;

         switch( vartype )
         {
            case SCIP_VARTYPE_BINARY :
               (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "b_%d", i);
               break;
            case SCIP_VARTYPE_INTEGER :
               (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "i_%d", i);
               break;
            case SCIP_VARTYPE_CONTINUOUS :
               (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "x_%d", i);
               break;
            default:
               SCIPABORT();
               break;
         }

         SCIP_CALL_THROW( SCIPcreateVarBasic(scip, &probdata->vars[i], name, -SCIPinfinity(scip), SCIPinfinity(scip), 0.0, vartype) );
         SCIP_CALL_THROW( SCIPaddVar(scip, probdata->vars[i]) );

         if( i < h.num_nl_vars_in_both + h.num_nl_vars_in_cons + h.num_nl_vars_in_objs )
         {
            SCIP_CALL_THROW( SCIPcreateExprVar(scip, &varexprs[i], probdata->vars[i], NULL, NULL) );
         }
      }

      // alloc some space for nonlinear constraints
      nlconslhss.resize(h.num_nl_cons, -SCIPinfinity(scip));
      nlconsrhss.resize(h.num_nl_cons,  SCIPinfinity(scip));
      nlconsexprs.resize(h.num_nl_cons, NULL);

      // create empty linear constraints
      linconss.reserve(h.num_algebraic_cons - h.num_nl_cons);
      for( int i = h.num_nl_cons; i < h.num_algebraic_cons; ++i )
      {
         (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "lc_%d", i);
         SCIP_CALL_THROW( SCIPcreateConsBasicLinear(scip, &cons, name, 0, NULL, NULL, -SCIPinfinity(scip), SCIPinfinity(scip)) );
         linconss.push_back(cons);
      }
   }

   SCIP_EXPR* OnNumber(
      double             value
      )
   {
      SCIP_EXPR* expr;

      SCIP_CALL_THROW( SCIPcreateExprValue(scip, &expr, value, NULL, NULL) );

      // remember that we have to release this expr
      exprstorelease.push_back(expr);

      return expr;
   }

   SCIP_EXPR* OnVariableRef(
      int                variableIndex
      )
   {
      assert(variableIndex >= 0);
      assert(variableIndex < (int)varexprs.size());
      assert(varexprs[variableIndex] != NULL);

      return varexprs[variableIndex];
   }

   SCIP_EXPR* OnUnary(
      mp::expr::Kind     kind,
      SCIP_EXPR*         child
      )
   {
      SCIP_EXPR* expr;
      switch( kind )
      {
         case mp::expr::MINUS:
         {
            SCIP_Real minusone = -1.0;
            SCIP_CALL_THROW( SCIPcreateExprSum(scip, &expr, 1, &child, &minusone, 0.0, NULL, NULL) );
            break;
         }

         case mp::expr::ABS:
            SCIP_CALL_THROW( SCIPcreateExprAbs(scip, &expr, child, NULL, NULL) );
            break;

         case mp::expr::POW2:
            SCIP_CALL_THROW( SCIPcreateExprPow(scip, &expr, child, 2.0, NULL, NULL) );
            break;

         case mp::expr::SQRT:
            SCIP_CALL_THROW( SCIPcreateExprPow(scip, &expr, child, 0.5, NULL, NULL) );
            break;

         case mp::expr::LOG:
            SCIP_CALL_THROW( SCIPcreateExprLog(scip, &expr, child, NULL, NULL) );
            break;

         case mp::expr::LOG10:  // 1/log(10)*log(child)
         {
            SCIP_EXPR* logexpr;
            SCIP_Real factor = 1.0/log(10.0);
            SCIP_CALL_THROW( SCIPcreateExprLog(scip, &logexpr, child, NULL, NULL) );
            SCIP_CALL_THROW( SCIPcreateExprSum(scip, &expr, 1, &child, &factor, 0.0, NULL, NULL) );
            SCIP_CALL_THROW( SCIPreleaseExpr(scip, &logexpr) );
            break;
         }

         case mp::expr::EXP:
            SCIP_CALL_THROW( SCIPcreateExprExp(scip, &expr, child, NULL, NULL) );
            break;

         case mp::expr::SIN:
            SCIP_CALL_THROW( SCIPcreateExprSin(scip, &expr, child, NULL, NULL) );
            break;

         case mp::expr::COS:
            SCIP_CALL_THROW( SCIPcreateExprCos(scip, &expr, child, NULL, NULL) );
            break;

         default:
            throw std::logic_error("Error: Unsupported AMPL function " + std::string(mp::expr::str(kind)));
      }

      // remember that we have to release this expr
      exprstorelease.push_back(expr);

      return expr;
   }

   SCIP_EXPR* OnBinary(
      mp::expr::Kind     kind,
      SCIP_EXPR*         firstChild,
      SCIP_EXPR*         secondChild
      )
   {
      SCIP_EXPR* expr;
      SCIP_EXPR* children[2] = { firstChild, secondChild };
      switch( kind )
      {
         case mp::expr::ADD:
            SCIP_CALL_THROW( SCIPcreateExprSum(scip, &expr, 2, children, NULL, 0.0, NULL, NULL) );
            break;

         case mp::expr::SUB:
         {
            SCIP_Real coefs[2] = { 1.0, -1.0 };
            SCIP_CALL_THROW( SCIPcreateExprSum(scip, &expr, 2, children, coefs, 0.0, NULL, NULL) );
            break;
         }

         case mp::expr::MUL:
            SCIP_CALL_THROW( SCIPcreateExprProduct(scip, &expr, 2, children, 1.0, NULL, NULL) );
            break;

         case mp::expr::DIV:
            SCIP_CALL_THROW( SCIPcreateExprPow(scip, &children[1], secondChild, -1.0, NULL, NULL) );
            SCIP_CALL_THROW( SCIPcreateExprProduct(scip, &expr, 2, children, 1.0, NULL, NULL) );
            SCIP_CALL_THROW( SCIPreleaseExpr(scip, &children[1]) );
            break;

         case mp::expr::POW:
         {
            // reformulate x^y as exp(y*log(x))
            SCIP_EXPR* prod;

            SCIP_CALL_THROW( SCIPcreateExprLog(scip, &children[0], firstChild, NULL, NULL) );  // log(firstChild)
            SCIP_CALL_THROW( SCIPcreateExprProduct(scip, &prod, 2, children, 1.0, NULL, NULL) );  // log(firstChild)*secondChild
            SCIP_CALL_THROW( SCIPcreateExprExp(scip, &expr, prod, NULL, NULL) );  // expr(log(firstChild)*secondChild)

            SCIP_CALL_THROW( SCIPreleaseExpr(scip, &prod) );
            SCIP_CALL_THROW( SCIPreleaseExpr(scip, &children[0]) );
            break;
         }

         case mp::expr::POW_CONST_BASE:
         {
            // reformulate constant^x as exp(y*log(constant)), require constant > 0.0
            SCIP_EXPR* prod;
            SCIP_Real coef;

            assert(SCIPisExprValue(scip, firstChild));
            if( SCIPgetValueExprValue(firstChild) <= 0.0 )
               throw std::logic_error("Error: Cannot handle pow(" + std::to_string(SCIPgetValueExprValue(firstChild)) + ", function)");

            coef = log(SCIPgetValueExprValue(firstChild)); // log(firstChild)
            SCIP_CALL_THROW( SCIPcreateExprSum(scip, &prod, 1, &secondChild, &coef, 0.0, NULL, NULL) );  // log(firstChild)*secondChild
            SCIP_CALL_THROW( SCIPcreateExprExp(scip, &expr, prod, NULL, NULL) );  // expr(log(firstChild)*secondChild)

            SCIP_CALL_THROW( SCIPreleaseExpr(scip, &prod) );
            break;
         }

         case mp::expr::POW_CONST_EXP:
            assert(SCIPisExprValue(scip, secondChild));
            SCIP_CALL_THROW( SCIPcreateExprPow(scip, &expr, firstChild, SCIPgetValueExprValue(secondChild), NULL, NULL) );
            break;

         default:
            throw std::logic_error("Error: Unsupported AMPL function " + std::string(mp::expr::str(kind)));
            break;
      }

      // remember that we have to release this expr
      exprstorelease.push_back(expr);

      return expr;
   }

   // Used for creating a list of terms in a sum
   // TODO copying around NumericArgHandler objects doesn't look efficient
   class NumericArgHandler : public std::vector<SCIP_EXPR*>
   {
   public:
      void AddArg(
         SCIP_EXPR*      term
         )
      {
         push_back(term);
      }
   };

   NumericArgHandler BeginSum(
      int                num_args
      )
   {
      NumericArgHandler h;
      h.reserve(num_args);
      return h;
   }

   SCIP_EXPR* EndSum(
      NumericArgHandler handler
      )
   {
      SCIP_EXPR* expr;
      SCIP_CALL_THROW( SCIPcreateExprSum(scip, &expr, (int)handler.size(), handler.data(), NULL, 0.0, NULL, NULL) );
      // remember that we have to release this expr
      exprstorelease.push_back(expr);
      return expr;
   }

   void OnObj(
      int                objectiveIndex,
      mp::obj::Type      type,
      SCIP_EXPR*         nonlinearExpression
      )
   {
      if( objectiveIndex >= 1 )
         throw std::logic_error("Error: Multiple objective functions not supported");

      SCIPsetObjsense(scip, type == mp::obj::Type::MAX ? SCIP_OBJSENSE_MAXIMIZE : SCIP_OBJSENSE_MINIMIZE);

      assert(objexpr == NULL);
      objexpr = nonlinearExpression;

      // FIXME how does AMPL handle a constant in a linear objective?
   }

   void OnAlgebraicCon(
      int                constraintIndex,
      SCIP_EXPR*         expr
      )
   {
      // nonlinear constraint iff expression is not NULL
      assert((constraintIndex < (int)nlconsexprs.size()) == (expr != NULL));
      if( expr != NULL )
         nlconsexprs[constraintIndex] = expr;
   }

   void OnVarBounds(
      int                variableIndex,
      double             variableLB,
      double             variableUB
      )
   {
      assert(variableIndex >= 0);
      assert(variableIndex < probdata->nvars);

      // as far as I see, ampl::mp gives -inf, +inf for no-bounds, which is always beyond SCIPinfinity()
      if( !SCIPisInfinity(scip, -variableLB) )
      {
         SCIP_CALL_THROW( SCIPchgVarLbGlobal(scip, probdata->vars[variableIndex], variableLB) );
      }
      if( !SCIPisInfinity(scip, variableUB) )
      {
         SCIP_CALL_THROW( SCIPchgVarUbGlobal(scip, probdata->vars[variableIndex], variableUB) );
      }
   }

   void OnConBounds(
      int                index,
      double             lb,
      double             ub
      )
   {
      assert(index >= 0);
      assert(index < (int)(nlconslhss.size() + linconss.size()));

      // nonlinear constraints are first
      if( index < (int)nlconslhss.size() )
      {
         if( !SCIPisInfinity(scip, -lb) )
            nlconslhss[index] = lb;
         if( !SCIPisInfinity(scip,  ub) )
            nlconsrhss[index] = ub;
      }
      else
      {
         if( !SCIPisInfinity(scip, -lb) )
         {
            SCIP_CALL_THROW( SCIPchgLhsLinear(scip, linconss[index - nlconslhss.size()], lb) );
         }
         if( !SCIPisInfinity(scip,  ub) )
         {
            SCIP_CALL_THROW( SCIPchgRhsLinear(scip, linconss[index - nlconslhss.size()], ub) );
         }
      }
   }

   class LinearPartHandler
   {
   private:
      SCIP* scip;
      SCIP_PROBDATA* probdata;
      SCIP_CONS* lincons;

   public:
      // constructor for linear constraint
      explicit LinearPartHandler(
         AMPLProblemHandler& amplph,
         int                 constraintIndex
         )
      : scip(amplph.scip),
        probdata(amplph.probdata),
        lincons(amplph.linconss[constraintIndex - amplph.nlconsexprs.size()])
      { }

      // constructor for linear objective
      explicit LinearPartHandler(
         AMPLProblemHandler& amplph
         )
      : scip(amplph.scip),
        probdata(amplph.probdata),
        lincons(NULL)
      { }

      void AddTerm(
         int             variableIndex,
         double          coefficient
         )
      {
         assert(variableIndex >= 0);
         assert(variableIndex < probdata->nvars);

         if( coefficient == 0.0 )
            return;

         if( lincons == NULL )
         {
            SCIP_CALL_THROW( SCIPchgVarObj(scip, probdata->vars[variableIndex], coefficient) );
         }
         else
         {
            SCIP_CALL_THROW( SCIPaddCoefLinear(scip, lincons, probdata->vars[variableIndex], coefficient) );
         }
      }
   };

   typedef LinearPartHandler LinearObjHandler;

   LinearPartHandler OnLinearObjExpr(
      int                objectiveIndex,
      int                /* numLinearTerms */
      )
   {
      if( objectiveIndex >= 1 )
         throw std::logic_error("Error: Multiple objective functions not supported");

      return LinearObjHandler(*this);
   }

   typedef LinearPartHandler LinearConHandler;

   LinearConHandler OnLinearConExpr(
      int                constraintIndex,
      int                /* numLinearTerms */
      )
   {
      return LinearConHandler(*this, constraintIndex);
   }

   void finalize()
   {
      char name[SCIP_MAXSTRLEN];

      // turn nonlinear objective into constraint
      // min f(x) -> min z s.t. f(x) - z <= 0
      // max f(x) -> max z s.t. 0 <= f(x) - z
      if( objexpr != NULL )
      {
         SCIP_CONS* objcons;
         SCIP_VAR* objvar;

         SCIP_CALL_THROW( SCIPcreateVarBasic(scip, &objvar, "objvar", -SCIPinfinity(scip), SCIPinfinity(scip), 1.0, SCIP_VARTYPE_CONTINUOUS) );
         SCIP_CALL_THROW( SCIPaddVar(scip, objvar) );

         SCIP_CALL_THROW( SCIPcreateConsBasicNonlinear(scip, &objcons, "objcons", objexpr,
            SCIPgetObjsense(scip) == SCIP_OBJSENSE_MINIMIZE ? -SCIPinfinity(scip) : 0.0,
            SCIPgetObjsense(scip) == SCIP_OBJSENSE_MAXIMIZE ?  SCIPinfinity(scip) : 0.0) );
         SCIP_CALL_THROW( SCIPaddLinearTermConsNonlinear(scip, objcons, -1.0, objvar) );
         SCIP_CALL_THROW( SCIPaddCons(scip, objcons) );

         SCIP_CALL_THROW( SCIPreleaseCons(scip, &objcons) );
         SCIP_CALL_THROW( SCIPreleaseVar(scip, &objvar) );
      }

      // create and add nonlinear constraints
      for( size_t i = 0; i < nlconsexprs.size(); ++i )
      {
         SCIP_CONS* cons;

         assert(nlconsexprs[i] != NULL);

         (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "nlc_%d", (int)i);
         SCIP_CALL_THROW( SCIPcreateConsBasicNonlinear(scip, &cons, name, nlconsexprs[i], nlconslhss[i], nlconsrhss[i]) );
         SCIP_CALL_THROW( SCIPaddCons(scip, cons) );
         SCIP_CALL_THROW( SCIPreleaseCons(scip, &cons) );
      }

      // add linear constraints
      for( size_t i = 0; i < linconss.size(); ++i )
      {
         SCIP_CALL_THROW( SCIPaddCons(scip, linconss[i]) );
         SCIP_CALL_THROW( SCIPreleaseCons(scip, &linconss[i]) );
      }

      // release created expressions (they should all be used in other expressions or constraints now)
      while( !exprstorelease.empty() )
      {
         SCIP_CALL_THROW( SCIPreleaseExpr(scip, &exprstorelease.back()) );
         exprstorelease.pop_back();
      }

      // release variable expressions (they should all be used in other expressions or constraints now)
      while( !varexprs.empty() )
      {
         SCIP_CALL_THROW( SCIPreleaseExpr(scip, &varexprs.back()) );
         varexprs.pop_back();
      }
   }
};


/*
 * Callback methods of probdata
 */

/** frees user data of original problem (called when the original problem is freed) */
static
SCIP_DECL_PROBDELORIG(probdataDelOrigNl)
{
   int i;

   assert((*probdata)->vars != NULL || (*probdata)->nvars == 0);

   for( i = 0; i < (*probdata)->nvars; ++i )
   {
      SCIP_CALL( SCIPreleaseVar(scip, &(*probdata)->vars[i]) );
   }
   SCIPfreeBlockMemoryArrayNull(scip, &(*probdata)->vars, (*probdata)->nvars);

   SCIPfreeMemory(scip, probdata);

   return SCIP_OKAY;
}

/*
 * Callback methods of reader
 */

/** copy method for reader plugins (called when SCIP copies plugins) */
static
SCIP_DECL_READERCOPY(readerCopyNl)
{  /*lint --e{715}*/
   assert(scip != NULL);

   SCIP_CALL( SCIPincludeReaderNl(scip) );

   return SCIP_OKAY;
}

/** problem reading method of reader */
static
SCIP_DECL_READERREAD(readerReadNl)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(reader != NULL);
   assert(filename != NULL);
   assert(result != NULL);

   try
   {
      AMPLProblemHandler handler(scip, filename);
      mp::ReadNLFile(filename, handler);
      handler.finalize();
   }
   catch( const std::exception& e )
   {
      // TODO distinguish exceptions and give different error return codes
      SCIPerrorMessage("Error when reading AMPL .nl file %s: %s", filename, e.what());
      return SCIP_READERROR;
   }

   *result = SCIP_SUCCESS;

   return SCIP_OKAY;
}

/*
 * reader specific interface methods
 */

/** includes the AMPL .nl file reader in SCIP */
SCIP_RETCODE SCIPincludeReaderNl(
   SCIP*                 scip                /**< SCIP data structure */
)
{
   SCIP_READER* reader = NULL;

   /* include reader */
   SCIP_CALL( SCIPincludeReaderBasic(scip, &reader, READER_NAME, READER_DESC, READER_EXTENSION, NULL) );
   assert(reader != NULL);

   /* set non fundamental callbacks via setter functions */
   SCIP_CALL( SCIPsetReaderCopy(scip, reader, readerCopyNl) );
   /* SCIP_CALL( SCIPsetReaderFree(scip, reader, readerFreeNl) ); */
   SCIP_CALL( SCIPsetReaderRead(scip, reader, readerReadNl) );

   SCIP_CALL( SCIPincludeExternalCodeInformation(scip, "AMPL/MP bb7d6166", "AMPL .nl file reader library (github.com/ampl/mp)") );

   return SCIP_OKAY;
}
