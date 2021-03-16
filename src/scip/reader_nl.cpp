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
 * For documentation on ampl::mp, see https://ampl.github.io and https://www.zverovich.net/2014/09/19/reading-nl-files.html.
 * For documentation on .nl files, see https://ampl.com/REFS/hooking2.pdf.
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

// a variant of SCIP_CALL that throws a std::logic_error if not SCIP_OKAY
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

/// problem data stored in SCIP
struct SCIP_ProbData
{
   SCIP_VAR**            vars;               /**< variables in the order given by AMPL */
   int                   nvars;              /**< number of variables */
};

/*
 * Local methods
 */

// forward declaration
static SCIP_DECL_PROBDELORIG(probdataDelOrigNl);

/// implementation of AMPL/MPs NLHandler that constructs a SCIP problem while a .nl file is read
class AMPLProblemHandler : public mp::NLHandler<AMPLProblemHandler, SCIP_EXPR*>
{
private:
   SCIP* scip;
   SCIP_PROBDATA* probdata;

   // variable expressions corresponding to nonlinear variables
   // created in OnHeader() and released in destructor
   // for reuse of var-expressions in OnVariableRef()
   std::vector<SCIP_EXPR*> varexprs;

   // constraints (collected here and added to SCIP in EndInput())
   // first nonliner, then linear
   std::vector<SCIP_CONS*> conss;

   // linear parts for nonlinear constraints
   // first collect and then add to constraints in EndInput()
   std::vector<std::vector<std::pair<SCIP_Real, SCIP_VAR*> > > nlconslin;

   // expression that represents a nonlinear objective function
   // used to create a corresponding constraint in EndInput(), unless NULL
   SCIP_EXPR* objexpr;

   // collect expressions that need to be released eventually
   // this are all expression that are returned to the AMPL/MP code in AMPLProblemHandler::OnXyz() functions
   // they need to be released exactly once, but after they are used in another expression or a constraint
   // as AMPL/MP may reuse expressions (common subexpressions), we don't release an expression when it is used
   // as a child or when constructing a constraint, but first collect them all and then release in destructor
   // alternatively, one could encapsulate SCIP_EXPR* into a small class that handles proper reference counting
   std::vector<SCIP_EXPR*> exprstorelease;

   // initial solution, if any
   SCIP_SOL* initsol;

   // opened files with column/variable and row/constraint names, or NULL
   fmt::File* colfile;
   fmt::File* rowfile;

   // get name from names strings, if possible
   // returns whether a name has been stored
   bool nextName(
      const char*&       namesbegin,         /**< current pointer into names string, or NULL */
      const char*        namesend,           /**< pointer to end of names string */
      char*              name                /**< buffer to store name, should have lenght SCIP_MAXSTRLEN */
   )
   {
      if( namesbegin == NULL )
         return false;

      // copy namesbegin into name until newline or namesend
      // updates namesbegin
      int nchars = 0;
      while( namesbegin != namesend )
      {
         if( nchars == SCIP_MAXSTRLEN )
         {
            SCIPverbMessage(scip, SCIP_VERBLEVEL_FULL, NULL, "name too long when parsing names file");
            // do no longer read names from this string (something seems awkward)
            namesbegin = NULL;
            return false;
         }
         if( *namesbegin == '\n' )
         {
            *name = '\0';
            ++namesbegin;
            return true;
         }
         *(name++) = *(namesbegin++);
         ++nchars;
      }

      SCIPverbMessage(scip, SCIP_VERBLEVEL_FULL, NULL, "missing newline when parsing names file");
      return false;
   }

public:
   /// constructor
   ///
   /// initializes SCIP problem and problem data
   AMPLProblemHandler(
      SCIP*              scip_,              /**< SCIP data structure */
      const char*        filename            /**< name of .nl file that is read */
      )
   : scip(scip_),
     probdata(NULL),
     objexpr(NULL),
     initsol(NULL),
     colfile(NULL),
     rowfile(NULL)
   {
      assert(scip != NULL);
      assert(filename != NULL);

      const char* filebasename = strrchr(const_cast<char*>(filename), '/');
      if( filebasename == NULL )
         filebasename = filename;
      else
         ++filebasename;

      SCIP_CALL_THROW( SCIPallocClearMemory(scip, &probdata) );

      // initialize empty SCIP problem
      SCIP_CALL_THROW( SCIPcreateProb(scip, filebasename, probdataDelOrigNl, NULL, NULL, NULL, NULL, NULL, probdata) );

      // open files with variable and constraint names
      int len = filename != NULL ? strlen(filename) : 0;
      if( len > 3 && filename[len-3] == '.' && filename[len-2] == 'n' && filename[len-1] == 'l' )
      {
         char* filename2;
         SCIP_CALL_THROW( SCIPallocBufferArray(scip, &filename2, len+1) );

         try
         {
            sprintf(filename2, "%.*s.col", len-3, filename);
            colfile = new fmt::File(filename2, fmt::File::RDONLY);

            filename2[len-2] = 'r';
            filename2[len-1] = 'o';
            filename2[len] = 'w';
            rowfile = new fmt::File(filename2, fmt::File::RDONLY);
         }
         catch( const fmt::SystemError& e )
         {
            // probably a file open error, probably because file not found
            // ignore, we can make up our own names
         }

         SCIPfreeBufferArray(scip, &filename2);
      }
   }

   /// destructor
   ///
   /// only asserts that cleanup() has been called, as we cannot throw an exception or return a SCIP_RETCODE here
   ~AMPLProblemHandler()
   {
      // exprs and linear constraint arrays should have been cleared up in cleanup()
      assert(conss.empty());
      assert(varexprs.empty());
      assert(exprstorelease.empty());

      delete colfile;
      delete rowfile;
   }

   /// process header of .nl files
   ///
   /// create and add variables, allocate constraints
   void OnHeader(
      const mp::NLHeader& h
      )
   {
      char name[SCIP_MAXSTRLEN];
      SCIP_CONS* cons;
      int nnlvars;

      assert(probdata->vars == NULL);

      // read variable and constraint names from file, if available, into memory
      // if not available, we will get varnamesbegin==NULL and consnamesbegin==NULL
      mp::MemoryMappedFile<> mapped_colfile;
      if( colfile != NULL )
         mapped_colfile.map(*colfile, "colfile");
      const char* varnamesbegin = mapped_colfile.start();
      const char* varnamesend = mapped_colfile.start() + mapped_colfile.size();

      mp::MemoryMappedFile<> mapped_rowfile;
      if( rowfile != NULL )
         mapped_rowfile.map(*rowfile, "rowfile");
      const char* consnamesbegin = mapped_rowfile.start();
      const char* consnamesend = mapped_rowfile.start() + mapped_rowfile.size();

      probdata->nvars = h.num_vars;
      SCIP_CALL_THROW( SCIPallocBlockMemoryArray(scip, &probdata->vars, probdata->nvars) );

      // number of nonlinear variables
      nnlvars = MAX(h.num_nl_vars_in_cons, h.num_nl_vars_in_objs);
      varexprs.resize(nnlvars);

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
         else if( i < h.num_nl_vars_in_cons - h.num_nl_integer_vars_in_cons )
            vartype = SCIP_VARTYPE_CONTINUOUS;
         else if( i < h.num_nl_vars_in_cons )
            vartype = SCIP_VARTYPE_INTEGER;
         // Nonlinear variables in objective
         else if( i < h.num_nl_vars_in_objs - h.num_nl_integer_vars_in_objs )
            vartype = SCIP_VARTYPE_CONTINUOUS;
         else if( i < h.num_nl_vars_in_objs )
            vartype = SCIP_VARTYPE_INTEGER;
         // Linear variables
         else if( i < h.num_vars - h.num_linear_binary_vars - h.num_linear_integer_vars )
            vartype = SCIP_VARTYPE_CONTINUOUS;
         else if( i < h.num_vars - h.num_linear_integer_vars )
            vartype = SCIP_VARTYPE_BINARY;
         else
            vartype = SCIP_VARTYPE_INTEGER;

         if( !nextName(varnamesbegin, varnamesend, name) )
         {
            // make up name if no names file or could not be read
            switch( vartype )
            {
               case SCIP_VARTYPE_BINARY :
                  (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "b%d", i);
                  break;
               case SCIP_VARTYPE_INTEGER :
                  (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "i%d", i);
                  break;
               case SCIP_VARTYPE_CONTINUOUS :
                  (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "x%d", i);
                  break;
               default:
                  SCIPABORT();
                  break;
            }
         }

         SCIP_CALL_THROW( SCIPcreateVarBasic(scip, &probdata->vars[i], name,
            vartype == SCIP_VARTYPE_BINARY ? 0.0 : -SCIPinfinity(scip),
            vartype == SCIP_VARTYPE_BINARY ? 1.0 :  SCIPinfinity(scip),
            0.0, vartype) );
         SCIP_CALL_THROW( SCIPaddVar(scip, probdata->vars[i]) );

         if( i < nnlvars )
         {
            SCIP_CALL_THROW( SCIPcreateExprVar(scip, &varexprs[i], probdata->vars[i], NULL, NULL) );
         }
      }

      // alloc some space for constraints
      conss.reserve(h.num_algebraic_cons);
      nlconslin.resize(h.num_nl_cons);

      // create empty nonlinear constraints
      // use expression == 0, because nonlinear constraint don't like to be without an expression
      SCIP_EXPR* dummyexpr;
      SCIP_CALL_THROW( SCIPcreateExprValue(scip, &dummyexpr, 0.0, NULL, NULL) );
      for( int i = 0; i < h.num_nl_cons; ++i )
      {
         // make up name if no names file or could not be read
         if( !nextName(consnamesbegin, consnamesend, name) )
            (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "nlc%d", i);

         SCIP_CALL_THROW( SCIPcreateConsBasicNonlinear(scip, &cons, name, dummyexpr, -SCIPinfinity(scip), SCIPinfinity(scip)) );
         conss.push_back(cons);
      }
      SCIP_CALL_THROW( SCIPreleaseExpr(scip, &dummyexpr) );

      // create empty linear constraints
      for( int i = h.num_nl_cons; i < h.num_algebraic_cons; ++i )
      {
         if( !nextName(consnamesbegin, consnamesend, name) )
            (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "lc%d", i);
         SCIP_CALL_THROW( SCIPcreateConsBasicLinear(scip, &cons, name, 0, NULL, NULL, -SCIPinfinity(scip), SCIPinfinity(scip)) );
         conss.push_back(cons);
      }
   }

   /// receive notification of a number in a nonlinear expression
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

   /// receive notification of a variable reference in a nonlinear expression
   SCIP_EXPR* OnVariableRef(
      int                variableIndex
      )
   {
      assert(variableIndex >= 0);
      assert(variableIndex < (int)varexprs.size());
      assert(varexprs[variableIndex] != NULL);

      return varexprs[variableIndex];
   }

   /// receive notification of a unary expression
   SCIP_EXPR* OnUnary(
      mp::expr::Kind     kind,
      SCIP_EXPR*         child
      )
   {
      SCIP_EXPR* expr;

      assert(child != NULL);

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
            OnUnhandled(mp::expr::str(kind));
            break;
      }

      // remember that we have to release this expr
      exprstorelease.push_back(expr);

      return expr;
   }

   /// receive notification of a binary expression
   SCIP_EXPR* OnBinary(
      mp::expr::Kind     kind,
      SCIP_EXPR*         firstChild,
      SCIP_EXPR*         secondChild
      )
   {
      SCIP_EXPR* expr;
      SCIP_EXPR* children[2] = { firstChild, secondChild };

      assert(firstChild != NULL);
      assert(secondChild != NULL);

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

         case mp::expr::POW_CONST_BASE:
         case mp::expr::POW_CONST_EXP:
         case mp::expr::POW:
            // with some .nl files, we seem to get mp::expr::POW even if base or exponent is constant,
            // so do not rely on kind but better check expr type
            if( SCIPisExprValue(scip, secondChild) )
            {
               SCIP_CALL_THROW( SCIPcreateExprPow(scip, &expr, firstChild, SCIPgetValueExprValue(secondChild), NULL, NULL) );
               break;
            }

            if( SCIPisExprValue(scip, firstChild) && SCIPgetValueExprValue(firstChild) > 0.0 )
            {
               // reformulate constant^x as exp(y*log(constant)), if constant > 0.0
               // if constant < 0, we create an expression and let cons_nonlinear figure out infeasibility somehow
               SCIP_EXPR* prod;

               SCIP_Real coef = log(SCIPgetValueExprValue(firstChild)); // log(firstChild)
               SCIP_CALL_THROW( SCIPcreateExprSum(scip, &prod, 1, &secondChild, &coef, 0.0, NULL, NULL) );  // log(firstChild)*secondChild
               SCIP_CALL_THROW( SCIPcreateExprExp(scip, &expr, prod, NULL, NULL) );  // expr(log(firstChild)*secondChild)

               SCIP_CALL_THROW( SCIPreleaseExpr(scip, &prod) );
               break;
            }

            {
               // reformulate x^y as exp(y*log(x))
               SCIP_EXPR* prod;

               assert(SCIPisExprValue(scip, secondChild));

               SCIP_CALL_THROW( SCIPcreateExprLog(scip, &children[0], firstChild, NULL, NULL) );  // log(firstChild)
               SCIP_CALL_THROW( SCIPcreateExprProduct(scip, &prod, 2, children, 1.0, NULL, NULL) );  // log(firstChild)*secondChild
               SCIP_CALL_THROW( SCIPcreateExprExp(scip, &expr, prod, NULL, NULL) );  // expr(log(firstChild)*secondChild)

               SCIP_CALL_THROW( SCIPreleaseExpr(scip, &prod) );
               SCIP_CALL_THROW( SCIPreleaseExpr(scip, &children[0]) );
               break;
            }

         default:
            OnUnhandled(mp::expr::str(kind));
            break;
      }

      // remember that we have to release this expr
      exprstorelease.push_back(expr);

      return expr;
   }

   /// handler to create a list of terms in a sum
   ///
   /// NumericArgHandler is copied around, so it keeps only a pointer (with reference counting) to actual data
   class NumericArgHandler
   {
   public:
      std::shared_ptr<std::vector<SCIP_EXPR*> > v;

      NumericArgHandler(
         int num_args
         )
      : v(new std::vector<SCIP_EXPR*>())
      {
         v->reserve(num_args);
      }

      void AddArg(
         SCIP_EXPR*      term
         )
      {
         v->push_back(term);
      }
   };

   /// receive notification of the beginning of a summation
   NumericArgHandler BeginSum(
      int                num_args
      )
   {
      NumericArgHandler h(num_args);
      return h;
   }

   /// receive notification of the end of a summation
   SCIP_EXPR* EndSum(
      NumericArgHandler handler
      )
   {
      SCIP_EXPR* expr;
      SCIP_CALL_THROW( SCIPcreateExprSum(scip, &expr, (int)handler.v->size(), handler.v->data(), NULL, 0.0, NULL, NULL) );
      // remember that we have to release this expr
      exprstorelease.push_back(expr);
      return expr;
   }

   /// receive notification of an objective type and the nonlinear part of an objective expression
   void OnObj(
      int                objectiveIndex,
      mp::obj::Type      type,
      SCIP_EXPR*         nonlinearExpression
      )
   {
      if( objectiveIndex >= 1 )
         OnUnhandled("multiple objective functions");

      SCIPsetObjsense(scip, type == mp::obj::Type::MAX ? SCIP_OBJSENSE_MAXIMIZE : SCIP_OBJSENSE_MINIMIZE);

      assert(objexpr == NULL);

      if( nonlinearExpression != NULL && SCIPisExprValue(scip, nonlinearExpression) )
      {
         // handle objective constant by adding a fixed variable for it
         SCIP_VAR* objconstvar;
         SCIP_Real objconst = SCIPgetValueExprValue(nonlinearExpression);

         SCIP_CALL_THROW( SCIPcreateVarBasic(scip, &objconstvar, "objconstant", objconst, objconst, 1.0, SCIP_VARTYPE_CONTINUOUS) );
         SCIP_CALL_THROW( SCIPaddVar(scip, objconstvar) );
         SCIP_CALL_THROW( SCIPreleaseVar(scip, &objconstvar) );
      }
      else
      {
         objexpr = nonlinearExpression;
      }
   }

   /// receive notification of an algebraic constraint expression
   void OnAlgebraicCon(
      int                constraintIndex,
      SCIP_EXPR*         expr
      )
   {
      if( expr != NULL )
      {
         SCIP_CALL_THROW( SCIPchgExprNonlinear(scip, conss[constraintIndex], expr) );
      }
   }

#if !1  // TODO
   LinearExprHandler BeginCommonExpr(int index, int num_linear_terms) {
     return builder_.common_expr(index).set_linear_expr(num_linear_terms);
   }
   void EndCommonExpr(int index, NumericExpr expr, int position) {
     SetCommonExpr(builder_.common_expr(index), expr, position);
   }

   Reference OnCommonExprRef(int expr_index) {
     return builder_.MakeCommonExpr(expr_index);
   }
#endif

   /// receive notification of variable bounds
   void OnVarBounds(
      int                variableIndex,
      double             variableLB,
      double             variableUB
      )
   {
      assert(variableIndex >= 0);
      assert(variableIndex < probdata->nvars);

      // as far as I see, ampl::mp gives -inf, +inf for no-bounds, which is always beyond SCIPinfinity()
      // we ignore bounds outside [-scipinfinity,scipinfinity] here
      // for binary variables, we also ignore bounds outside [0,1]
      if( variableLB > (SCIPvarGetType(probdata->vars[variableIndex]) == SCIP_VARTYPE_BINARY ? 0.0 : -SCIPinfinity(scip)) )
      {
         SCIP_CALL_THROW( SCIPchgVarLbGlobal(scip, probdata->vars[variableIndex], variableLB) );
      }
      if( variableUB < (SCIPvarGetType(probdata->vars[variableIndex]) == SCIP_VARTYPE_BINARY ? 1.0 :  SCIPinfinity(scip)) )
      {
         SCIP_CALL_THROW( SCIPchgVarUbGlobal(scip, probdata->vars[variableIndex], variableUB) );
      }
   }

   /// receive notification of constraint sides
   void OnConBounds(
      int                index,
      double             lb,
      double             ub
      )
   {
      assert(index >= 0);
      assert(index < (int)conss.size());

      // nonlinear constraints are first
      if( index < (int)nlconslin.size() )
      {
         if( !SCIPisInfinity(scip, -lb) )
         {
            SCIP_CALL_THROW( SCIPchgLhsNonlinear(scip, conss[index], lb) );
         }
         if( !SCIPisInfinity(scip,  ub) )
         {
            SCIP_CALL_THROW( SCIPchgRhsNonlinear(scip, conss[index], ub) );
         }
      }
      else
      {
         if( !SCIPisInfinity(scip, -lb) )
         {
            SCIP_CALL_THROW( SCIPchgLhsLinear(scip, conss[index], lb) );
         }
         if( !SCIPisInfinity(scip,  ub) )
         {
            SCIP_CALL_THROW( SCIPchgRhsLinear(scip, conss[index], ub) );
         }
      }
   }

   /// receive notification of the initial value for a variable
   void OnInitialValue(
      int                var_index,
      double             value
      )
   {
      if( initsol == NULL )
      {
         SCIP_CALL_THROW( SCIPcreateSol(scip, &initsol, NULL) );
      }

      SCIP_CALL_THROW( SCIPsetSolVal(scip, initsol, probdata->vars[var_index], value) );
   }

   /// receives notification of the initial value for a dual variable
   void OnInitialDualValue(
      int                /* con_index */,
      double             /* value */
      )
   {
      // ignore initial dual value
   }

   /// receives notification of Jacobian column sizes
   ColumnSizeHandler OnColumnSizes()
   {
      /// use ColumnSizeHandler from upper class, which does nothing
      return ColumnSizeHandler();
   }

   // TODO finish
   template<typename T> class SuffixHandler
   {
   private:
      AMPLProblemHandler& amplph;
   public:
      SuffixHandler(
         AMPLProblemHandler& amplph_,
         fmt::StringRef      name,
         mp::suf::Kind       kind,
         int                 num_values
         )
      : amplph(amplph_)
      {
         switch( kind )
         {
            case mp::suf::Kind::CON:
               if( strncmp(name.data(), "initial", name.size()) == 0 )
               {

               }
               else if( strncmp(name.data(), "separate", name.size()) == 0 )
               {

               }
               else if( strncmp(name.data(), "enforce", name.size()) == 0 )
               {

               }
               else if( strncmp(name.data(), "check", name.size()) == 0 )
               {

               }
               else if( strncmp(name.data(), "propagate", name.size()) == 0 )
               {

               }
               else if( strncmp(name.data(), "dynamic", name.size()) == 0 )
               {

               }
               else if( strncmp(name.data(), "removable", name.size()) == 0 )
               {

               }
               else if( strncmp(name.data(), "sos", name.size()) == 0 )
               {
                  // SOS
               }
               SCIPverbMessage(amplph.scip, SCIP_VERBLEVEL_HIGH, NULL, "Unknown constraint suffix <%.*s>. Ignoring.\n", name.size(), name.data());
               break;

            case mp::suf::Kind::VAR:
            {
               if( strncmp(name.data(), "initial", name.size()) == 0 )
               {

               }
               else if( strncmp(name.data(), "removable", name.size()) == 0 )
               {

               }
               else if( strncmp(name.data(), "ref", name.size()) == 0 )
               {
                  // SOS, real
               }
               else if( strncmp(name.data(), "sos", name.size()) == 0 )
               {
                  // SOS
               }
               else if( strncmp(name.data(), "sosno", name.size()) == 0 )
               {
                  // SOS, real (?)
               }
               else if( strncmp(name.data(), "priority", name.size()) == 0 )
               {
                  // SOS, not real (?)
               }
               SCIPverbMessage(amplph.scip, SCIP_VERBLEVEL_HIGH, NULL, "Unknown variable suffix <%.*s>. Ignoring.\n", name.size(), name.data());
               break;

            case mp::suf::Kind::OBJ:
               SCIPverbMessage(amplph.scip, SCIP_VERBLEVEL_HIGH, NULL, "Unknown objective suffix <%.*s>. Ignoring.\n", name.size(), name.data());
               break;

            case mp::suf::Kind::PROBLEM:
               SCIPverbMessage(amplph.scip, SCIP_VERBLEVEL_HIGH, NULL, "Unknown problem suffix <%.*s>. Ignoring.\n", name.size(), name.data());
               break;
            }
         }
      }

      void SetValue(
         int index,
         T   value
      )
      { }
   };

   typedef SuffixHandler<int> IntSuffixHandler;
   /// receive notification of an integer suffix
   IntSuffixHandler OnIntSuffix(
      fmt::StringRef     name,               /**< suffix name, not null-terminated */
      mp::suf::Kind      kind,               /**< suffix kind */
      int                num_values
      )
   {
      return IntSuffixHandler(*this, name, kind, num_values);
   }

   typedef SuffixHandler<SCIP_Real> DblSuffixHandler;
   /// receive notification of a double suffix
   DblSuffixHandler OnDblSuffix(
      fmt::StringRef     name,               /**< suffix name, not null-terminated */
      mp::suf::Kind      kind,               /**< suffix kind */
      int                num_values
      )
   {
      return DblSuffixHandler(*this, name, kind, num_values);
   }

   /// handles receiving the linear part of an objective or constraint
   ///
   /// for objective, set the objective-coefficient of the variable
   /// for linear constraints, add to the constraint
   /// for nonlinear constraints, add to nlconslin vector; adding to constraint later
   class LinearPartHandler
   {
   private:
      AMPLProblemHandler& amplph;
      int constraintIndex;

   public:
      // constructor for constraint
      explicit LinearPartHandler(
         AMPLProblemHandler& amplph_,
         int                 constraintIndex_
         )
      : amplph(amplph_),
        constraintIndex(constraintIndex_)
      {
         assert(constraintIndex_ >= 0);
         assert(constraintIndex_ < (int)amplph.conss.size());
      }

      // constructor for linear objective
      explicit LinearPartHandler(
         AMPLProblemHandler& amplph_
         )
      : amplph(amplph_),
        constraintIndex(-1)
      { }

      void AddTerm(
         int             variableIndex,
         double          coefficient
         )
      {
         assert(variableIndex >= 0);
         assert(variableIndex < amplph.probdata->nvars);

         if( coefficient == 0.0 )
            return;

         if( constraintIndex < 0 )
         {
            SCIP_CALL_THROW( SCIPchgVarObj(amplph.scip, amplph.probdata->vars[variableIndex], coefficient) );
         }
         else if( constraintIndex < (int)amplph.nlconslin.size() )
         {
            amplph.nlconslin[constraintIndex].push_back(std::pair<SCIP_Real, SCIP_VAR*>(coefficient, amplph.probdata->vars[variableIndex]));
         }
         else
         {
            SCIP_CONS* lincons = amplph.conss.at(constraintIndex);
            SCIP_CALL_THROW( SCIPaddCoefLinear(amplph.scip, lincons, amplph.probdata->vars[variableIndex], coefficient) );
         }
      }
   };

   typedef LinearPartHandler LinearObjHandler;

   /// receive notification of the linear part of an objective
   LinearPartHandler OnLinearObjExpr(
      int                objectiveIndex,
      int                /* numLinearTerms */
      )
   {
      if( objectiveIndex >= 1 )
         OnUnhandled("multiple objective functions");

      return LinearObjHandler(*this);
   }

   typedef LinearPartHandler LinearConHandler;

   /// receive notification of the linear part of a constraint
   LinearConHandler OnLinearConExpr(
      int                constraintIndex,
      int                /* numLinearTerms */
      )
   {
      return LinearConHandler(*this, constraintIndex);
   }

   /// receive notification of the end of the input
   ///
   /// - setup all nonlinear constraints and add them to SCIP
   /// - add linear constraints to SCIP (should be after nonlinear ones to respect order in .nl file)
   /// - add initial solution, if initial values were given
   void EndInput()
   {
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

      // add linear terms to expressions of nonlinear constraints (should be ok to do this one-by-one for now)
      for( size_t i = 0; i < nlconslin.size(); ++i )
      {
         for( size_t j = 0; j < nlconslin[i].size(); ++j )
         {
            SCIP_CALL_THROW( SCIPaddLinearTermConsNonlinear(scip, conss[i], nlconslin[i][j].first, nlconslin[i][j].second) );
         }
      }

      // add and release constraints
      for( size_t i = 0; i < conss.size(); ++i )
      {
         SCIP_CALL_THROW( SCIPaddCons(scip, conss[i]) );
         SCIP_CALL_THROW( SCIPreleaseCons(scip, &conss[i]) );
      }
      conss.clear();

      // add initial solution
      if( initsol != NULL )
      {
         SCIP_Bool stored;
         SCIP_CALL_THROW( SCIPaddSolFree(scip, &initsol, &stored) );
      }

      // release expressions
      SCIP_CALL_THROW( cleanup() );
   }

   /// releases expressions and linear constraints from data
   ///
   /// should be called if there was an error while reading the .nl file
   /// this is not in the destructor, because we want to return SCIP_RETCODE
   SCIP_RETCODE cleanup()
   {
      // release constraints (in case EndInput() wasn't called)
      while( !conss.empty() )
      {
         SCIP_CALL( SCIPreleaseCons(scip, &conss.back()) );
         conss.pop_back();
      }

      // release initial sol (in case EndInput() wasn't called)
      if( initsol != NULL )
      {
         SCIP_CALL( SCIPfreeSol(scip, &initsol) );
      }

      // release created expressions (they should all be used in other expressions or constraints now)
      while( !exprstorelease.empty() )
      {
         SCIP_CALL( SCIPreleaseExpr(scip, &exprstorelease.back()) );
         exprstorelease.pop_back();
      }

      // release variable expressions (they should all be used in other expressions or constraints now)
      while( !varexprs.empty() )
      {
         SCIP_CALL( SCIPreleaseExpr(scip, &varexprs.back()) );
         varexprs.pop_back();
      }

      return SCIP_OKAY;
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
      // try to read the .nl file and setup SCIP problem
      AMPLProblemHandler handler(scip, filename);
      try
      {
         mp::ReadNLFile(filename, handler);
      }
      catch( const mp::UnsupportedError& e )
      {
         SCIPerrorMessage("unsupported construct in AMPL .nl file %s: %s\n", filename, e.what());

         SCIP_CALL( handler.cleanup() );

         return SCIP_READERROR;
      }
      catch( const fmt::SystemError& e )
      {
         // probably a file open error, probably because file not found
         SCIPerrorMessage("%s\n", e.what());

         SCIP_CALL( handler.cleanup() );

         return SCIP_NOFILE;
      }
      catch( const std::bad_alloc& e )
      {
         SCIPerrorMessage("Out of memory: %s\n", e.what());

         SCIP_CALL( handler.cleanup() );

         return SCIP_NOMEMORY;
      }
   }
   catch( const std::exception& e )
   {
      SCIPerrorMessage("%s\n", e.what());
      return SCIP_ERROR;
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
