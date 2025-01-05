/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*  Copyright (c) 2002-2025 Zuse Institute Berlin (ZIB)                      */
/*                                                                           */
/*  Licensed under the Apache License, Version 2.0 (the "License");          */
/*  you may not use this file except in compliance with the License.         */
/*  You may obtain a copy of the License at                                  */
/*                                                                           */
/*      http://www.apache.org/licenses/LICENSE-2.0                           */
/*                                                                           */
/*  Unless required by applicable law or agreed to in writing, software      */
/*  distributed under the License is distributed on an "AS IS" BASIS,        */
/*  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. */
/*  See the License for the specific language governing permissions and      */
/*  limitations under the License.                                           */
/*                                                                           */
/*  You should have received a copy of the Apache-2.0 license                */
/*  along with SCIP; see the file LICENSE. If not visit scipopt.org.         */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   reader_nl.cpp
 * @ingroup DEFPLUGINS_READER
 * @brief  AMPL .nl file reader
 * @author Stefan Vigerske
 *
 * For documentation on ampl::mp, see https://ampl.github.io and https://www.zverovich.net/2014/09/19/reading-nl-files.html.
 * For documentation on .nl files, see https://ampl.com/REFS/hooking2.pdf.
 */

/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <string>
#include <vector>
#include <map>

#include "scip/reader_nl.h"
#include "scip/cons_linear.h"
#include "scip/cons_setppc.h"
#include "scip/cons_logicor.h"
#include "scip/cons_knapsack.h"
#include "scip/cons_varbound.h"
#include "scip/cons_nonlinear.h"
#include "scip/cons_sos1.h"
#include "scip/cons_sos2.h"
#include "scip/cons_and.h"
#include "scip/cons_or.h"
#include "scip/cons_xor.h"
#include "scip/expr_var.h"
#include "scip/expr_value.h"
#include "scip/expr_sum.h"
#include "scip/expr_product.h"
#include "scip/expr_pow.h"
#include "scip/expr_log.h"
#include "scip/expr_exp.h"
#include "scip/expr_trig.h"
#include "scip/expr_abs.h"

// disable -Wshadow warnings for upcoming includes of AMPL/MP
// disable -Wimplicit-fallthrough as I don't want to maintain extra comments in AMPL/MP code to suppress these
#ifdef __GNUC__
#pragma GCC diagnostic ignored "-Wshadow"
#if __GNUC__ >= 7
#pragma GCC diagnostic ignored "-Wimplicit-fallthrough"
#endif
#endif

#include "mp/nl-reader.h"
#include "mp/nl-writer2.hpp"
#include "mp/nl-opcodes.h"

#define READER_NAME             "nlreader"
#define READER_DESC             "AMPL .nl file reader"
#define READER_EXTENSION        "nl"

// a variant of SCIP_CALL that throws a std::logic_error if not SCIP_OKAY
// (using cast to long long to work around issues with old MSVC)
#define SCIP_CALL_THROW(x) \
   do                                                                                                   \
   {                                                                                                    \
      SCIP_RETCODE throw_retcode;                                                                       \
      if( ((throw_retcode) = (x)) != SCIP_OKAY )                                                        \
         throw std::logic_error("Error <" + std::to_string((long long)throw_retcode) + "> in function call at reader_nl.cpp:" + std::to_string(__LINE__)); \
   }                                                                                                    \
   while( false )

/*
 * Data structures
 */

/// problem data stored in SCIP
struct SCIP_ProbData
{
   char*                 filenamestub;       /**< name of input file, without .nl extension; array is long enough to hold 5 extra chars */
   int                   filenamestublen;    /**< length of filenamestub string */

   int                   amplopts[mp::MAX_AMPL_OPTIONS];  /**< AMPL options from .nl header */
   int                   namplopts;          /**< number of AMPL options from .nl header */

   SCIP_VAR**            vars;               /**< variables in the order given by AMPL */
   int                   nvars;              /**< number of variables */

   SCIP_CONS**           conss;              /**< constraints in the order given by AMPL */
   int                   nconss;             /**< number of constraints */

   SCIP_Bool             islp;               /**< whether problem is an LP (only linear constraints, only continuous vars) */
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

   // linear parts for nonlinear constraints
   // first collect and then add to constraints in EndInput()
   std::vector<std::vector<std::pair<SCIP_Real, SCIP_VAR*> > > nlconslin;

   // expression that represents a nonlinear objective function
   // used to create a corresponding constraint in EndInput(), unless NULL
   SCIP_EXPR* objexpr;

   // common expressions (defined variables from statements like "var xsqr = x^2;" in an AMPL model)
   // they are constructed by BeginCommonExpr/EndCommonExpr below and are referenced by index in OnCommonExprRef
   std::vector<SCIP_EXPR*> commonexprs;

   // collect expressions that need to be released eventually
   // this are all expression that are returned to the AMPL/MP code in AMPLProblemHandler::OnXyz() functions
   // they need to be released exactly once, but after they are used in another expression or a constraint
   // as AMPL/MP may reuse expressions (common subexpressions), we don't release an expression when it is used
   // as a child or when constructing a constraint, but first collect them all and then release in destructor
   // alternatively, one could encapsulate SCIP_EXPR* into a small class that handles proper reference counting
   std::vector<SCIP_EXPR*> exprstorelease;

   // count on variables or constraints added for logical expressions
   int logiccount;

   // SOS constraints
   // collected while handling suffixes in SuffixHandler
   // sosvars maps the SOS index (can be negative) to the indices of the variables in the SOS
   // sosweights gives for each variable its weight in the SOS it appears in (if any)
   std::map<int, std::vector<int> > sosvars;
   std::vector<int> sosweights;

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
      char*              name                /**< buffer to store name, should have length SCIP_MAXSTRLEN */
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

   /// returns variable or value for given expression
   ///
   /// if expression is variable, ensure that it is a binary variable and set var
   /// if expression is value, then set val to whether value is nonzero and set var to NULL
   /// otherwise throw UnsupportedError exception
   void LogicalExprToVarVal(
      LogicalExpr        expr,
      SCIP_VAR*&         var,
      SCIP_Bool&         val
      )
   {
      assert(expr != NULL);

      if( SCIPisExprVar(scip, expr) )
      {
         var = SCIPgetVarExprVar(expr);
         if( SCIPvarGetType(var) != SCIP_VARTYPE_BINARY )
         {
            SCIP_Bool infeas;
            SCIP_Bool tightened;
            SCIP_CALL_THROW( SCIPchgVarType(scip, var, SCIP_VARTYPE_BINARY, &infeas) );
            assert(!infeas);
            SCIP_CALL_THROW( SCIPtightenVarLbGlobal(scip, var, 0.0, TRUE, &infeas, &tightened) );
            assert(!infeas);
            SCIP_CALL_THROW( SCIPtightenVarUbGlobal(scip, var, 1.0, TRUE, &infeas, &tightened) );
            assert(!infeas);
         }
         val = FALSE;  // for scan-build

         return;
      }

      if( SCIPisExprValue(scip, expr) )
      {
         var = NULL;
         val = SCIPgetValueExprValue(expr) != 0.0;
         return;
      }

      OnUnhandled("logical expression must be binary or constant");
   }

public:
   /// constructor
   ///
   /// initializes SCIP problem and problem data
   AMPLProblemHandler(
      SCIP*              scip_,              ///< SCIP data structure
      const char*        filename            ///< name of .nl file that is read
      )
   : scip(scip_),
     probdata(NULL),
     objexpr(NULL),
     logiccount(0),
     initsol(NULL),
     colfile(NULL),
     rowfile(NULL)
   {
      assert(scip != NULL);
      assert(filename != NULL);

      SCIP_CALL_THROW( SCIPallocClearMemory(scip, &probdata) );

      /* get name of input file without file extension (if any) */
      const char* extstart = strrchr(const_cast<char*>(filename), '.');
      if( extstart != NULL )
         probdata->filenamestublen = extstart - filename;
      else
         probdata->filenamestublen = strlen(filename);
      assert(probdata->filenamestublen > 0);
      SCIP_CALL_THROW( SCIPallocBlockMemoryArray(scip, &probdata->filenamestub, probdata->filenamestublen + 5) );
      memcpy(probdata->filenamestub, filename, probdata->filenamestublen);
      probdata->filenamestub[probdata->filenamestublen] = '\0';

      /* derive probname from name of input file without path and extension */
      const char* probname = strrchr(probdata->filenamestub, '/');
      if( probname == NULL )
         probname = probdata->filenamestub;
      else
         ++probname;

      // initialize empty SCIP problem
      SCIP_CALL_THROW( SCIPcreateProb(scip, probname, probdataDelOrigNl, NULL, NULL, NULL, NULL, NULL, probdata) );

      // try to open files with variable and constraint names
      // temporarily add ".col" and ".row", respectively, to filenamestub
      try
      {
         probdata->filenamestub[probdata->filenamestublen] = '.';
         probdata->filenamestub[probdata->filenamestublen+1] = 'c';
         probdata->filenamestub[probdata->filenamestublen+2] = 'o';
         probdata->filenamestub[probdata->filenamestublen+3] = 'l';
         probdata->filenamestub[probdata->filenamestublen+4] = '\0';
         colfile = new fmt::File(probdata->filenamestub, fmt::File::RDONLY);

         probdata->filenamestub[probdata->filenamestublen+1] = 'r';
         probdata->filenamestub[probdata->filenamestublen+3] = 'w';
         rowfile = new fmt::File(probdata->filenamestub, fmt::File::RDONLY);
      }
      catch( const fmt::SystemError& e )
      {
         // probably a file open error, probably because file not found
         // ignore, we can make up our own names
      }
      probdata->filenamestub[probdata->filenamestublen] = '\0';
   }

   AMPLProblemHandler(const AMPLProblemHandler&) = delete;
   AMPLProblemHandler& operator=(const AMPLProblemHandler&) = delete;

   /// destructor
   ///
   /// only asserts that cleanup() has been called, as we cannot throw an exception or return a SCIP_RETCODE here
   ~AMPLProblemHandler()
   {
      // exprs and linear constraint arrays should have been cleared up in cleanup()
      assert(varexprs.empty());
      assert(exprstorelease.empty());

      delete colfile;
      delete rowfile;
   }

   /// process header of .nl files
   ///
   /// create and add variables, allocate constraints
   void OnHeader(
      const mp::NLHeader& h                  ///< header data
      )
   {
      char name[SCIP_MAXSTRLEN];
      int nnlvars;

      assert(probdata->vars == NULL);
      assert(probdata->conss == NULL);

      probdata->namplopts = h.num_ampl_options;
      BMScopyMemoryArray(probdata->amplopts, h.ampl_options, h.num_ampl_options);

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
               // coverity[deadcode]
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

      // alloc some space for algebraic constraints
      probdata->nconss = h.num_algebraic_cons;
      SCIP_CALL_THROW( SCIPallocBlockMemoryArray(scip, &probdata->conss, probdata->nconss) );
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

         SCIP_CALL_THROW( SCIPcreateConsBasicNonlinear(scip, &probdata->conss[i], name, dummyexpr, -SCIPinfinity(scip), SCIPinfinity(scip)) );
      }
      SCIP_CALL_THROW( SCIPreleaseExpr(scip, &dummyexpr) );

      // create empty linear constraints
      for( int i = h.num_nl_cons; i < h.num_algebraic_cons; ++i )
      {
         if( !nextName(consnamesbegin, consnamesend, name) )
            (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "lc%d", i);
         SCIP_CALL_THROW( SCIPcreateConsBasicLinear(scip, &probdata->conss[i], name, 0, NULL, NULL, -SCIPinfinity(scip), SCIPinfinity(scip)) );
      }

      if( h.num_nl_cons == 0 && h.num_logical_cons == 0 && h.num_integer_vars() == 0 )
         probdata->islp = true;

      // alloc space for common expressions
      commonexprs.resize(h.num_common_exprs());
   }

   /// receive notification of a number in a nonlinear expression
   SCIP_EXPR* OnNumber(
      double             value               ///< value
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
      int                variableIndex       ///< AMPL index of variable
      )
   {
      assert(variableIndex >= 0);
      assert(variableIndex < (int)varexprs.size());
      assert(varexprs[variableIndex] != NULL);

      return varexprs[variableIndex];
   }

   /// receive notification of a unary expression
   SCIP_EXPR* OnUnary(
      mp::expr::Kind     kind,               ///< expression operator
      SCIP_EXPR*         child               ///< argument
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
            SCIP_CALL_THROW( SCIPcreateExprSum(scip, &expr, 1, &logexpr, &factor, 0.0, NULL, NULL) );
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
            return NULL;
      }

      // remember that we have to release this expr
      exprstorelease.push_back(expr);

      return expr;
   }

   /// receive notification of a binary expression
   SCIP_EXPR* OnBinary(
      mp::expr::Kind     kind,               ///< expression operand
      SCIP_EXPR*         firstChild,         ///< first argument
      SCIP_EXPR*         secondChild         ///< second argument
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
               // reformulate constant^y as exp(y*log(constant)), if constant > 0.0
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
            return NULL;
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

      /// constructor
      NumericArgHandler(
         int             num_args            ///< number of terms to expect
         )
      : v(new std::vector<SCIP_EXPR*>())
      {
         v->reserve(num_args);
      }

      /// adds term to sum
      void AddArg(
         SCIP_EXPR*      term                ///< term to add
         )
      {
         v->push_back(term);
      }
   };

   /// receive notification of the beginning of a summation
   NumericArgHandler BeginSum(
      int                num_args            ///< number of terms to expect
      )
   {
      NumericArgHandler h(num_args);
      return h;
   }

   /// receive notification of the end of a summation
   SCIP_EXPR* EndSum(
      NumericArgHandler  handler             ///< handler that handled the sum
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
      int                objectiveIndex,     ///< index of objective
      mp::obj::Type      type,               ///< objective sense
      SCIP_EXPR*         nonlinearExpression ///< nonlinear part of objective function
      )
   {
      if( objectiveIndex >= 1 )
         OnUnhandled("multiple objective functions");

      SCIP_CALL_THROW( SCIPsetObjsense(scip, type == mp::obj::Type::MAX ? SCIP_OBJSENSE_MAXIMIZE : SCIP_OBJSENSE_MINIMIZE) );

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
      int                constraintIndex,    ///< index of constraint
      SCIP_EXPR*         expr                ///< nonlinear part of constraint
      )
   {
      if( expr != NULL )
      {
         SCIP_CALL_THROW( SCIPchgExprNonlinear(scip, probdata->conss[constraintIndex], expr) );
      }
   }

   /// receives notification of a logical constraint expression
   void OnLogicalCon(
      int                index,
      LogicalExpr        expr
      )
   {
      if( expr != NULL )
      {
         SCIP_CONS* cons;
         SCIP_CALL_THROW( SCIPcreateConsBasicNonlinear(scip, &cons, "logiccons", expr, 1.0, 1.0) );
         SCIP_CALL_THROW( SCIPaddCons(scip, cons) );
         SCIP_CALL_THROW( SCIPreleaseCons(scip, &cons) );
      }
   }

   /// handles linear part of a common expression
   /// sets up a sum expression, if the linear part isn't empty
   class LinearExprHandler
   {
   private:
      AMPLProblemHandler& amplph;
      SCIP_EXPR*          commonexpr;

   public:
      /// constructor
      LinearExprHandler(
         AMPLProblemHandler& amplph_,        ///< problem handler
         int                 index,          ///< index of common expression
         int                 num_linear_terms///< number of terms to expect
         )
      : amplph(amplph_),
        commonexpr(NULL)
      {
         if( num_linear_terms > 0 )
         {
            SCIP_CALL_THROW( SCIPcreateExprSum(amplph.scip, &commonexpr, 0, NULL, NULL, 0.0, NULL, NULL) );
            amplph.commonexprs[index] = commonexpr;
            amplph.exprstorelease.push_back(commonexpr);
         }
      }

      /// receives notification of a term in the linear expression
      void AddTerm(
         int             var_index,          ///< AMPL index of variable
         double          coef                ///< variable coefficient
         )
      {
         assert(commonexpr != NULL);

         if( coef == 0.0 )
            return;

         if( var_index < (int)amplph.varexprs.size() )
         {
            SCIP_CALL_THROW( SCIPappendExprSumExpr(amplph.scip, commonexpr, amplph.varexprs[var_index], coef) );
         }
         else
         {
            // the index variable is linear (not sure this can happen here)
            assert(var_index < amplph.probdata->nvars);
            SCIP_EXPR* varexpr;
            SCIP_CALL_THROW( SCIPcreateExprVar(amplph.scip, &varexpr, amplph.probdata->vars[var_index], NULL, NULL) );
            SCIP_CALL_THROW( SCIPappendExprSumExpr(amplph.scip, commonexpr, varexpr, coef) );
            SCIP_CALL_THROW( SCIPreleaseExpr(amplph.scip, &varexpr) );
         }
      }
   };

   /// receive notification of the beginning of a common expression (defined variable)
   LinearExprHandler BeginCommonExpr(
      int                index,              ///< index of common expression
      int                num_linear_terms    ///< number of terms to expect
      )
   {
      assert(index >= 0);
      assert(index < (int)commonexprs.size());

      return LinearExprHandler(*this, index, num_linear_terms);
   }

   /// receive notification of the end of a common expression
   void EndCommonExpr(
      int                index,              ///< index of common expression
      SCIP_EXPR*         expr,               ///< nonlinear part of common expression
      int                /* position */      ///< argument that doesn't seem to have any purpose
      )
   {
      if( commonexprs[index] != NULL )
      {
         // add expr, if any, to linear part
         if( expr != NULL )
         {
            SCIP_CALL_THROW( SCIPappendExprSumExpr(scip, commonexprs[index], expr, 1.0) );
         }
      }
      else if( expr != NULL )
      {
         commonexprs[index] = expr;
      }
   }

   /// receive notification of a common expression (defined variable) reference
   SCIP_EXPR* OnCommonExprRef(
      int                expr_index          ///< index of common expression
      )
   {
      assert(expr_index >= 0);
      assert(expr_index < (int)commonexprs.size());
      assert(commonexprs[expr_index] != NULL);
      return commonexprs[expr_index];
   }

   /// receive notification of variable bounds
   void OnVarBounds(
      int                variableIndex,      ///< AMPL index of variable
      double             variableLB,         ///< variable lower bound
      double             variableUB          ///< variable upper bound
      )
   {
      assert(variableIndex >= 0);
      assert(variableIndex < probdata->nvars);

      // as far as I see, ampl::mp gives -inf, +inf for no-bounds, which is always beyond SCIPinfinity()
      // we ignore bounds outside [-scipinfinity,scipinfinity] here
      // for binary variables, we also ignore bounds outside [0,1]
      SCIP_Bool binary = (SCIPvarGetType(probdata->vars[variableIndex]) == SCIP_VARTYPE_BINARY);
      if( variableLB > (binary ? 0.0 : -SCIPinfinity(scip)) )
      {
         SCIP_CALL_THROW( SCIPchgVarLbGlobal(scip, probdata->vars[variableIndex], variableLB) );
      }
      if( variableUB < (binary ? 1.0 :  SCIPinfinity(scip)) )
      {
         SCIP_CALL_THROW( SCIPchgVarUbGlobal(scip, probdata->vars[variableIndex], variableUB) );
      }
   }

   /// receive notification of constraint sides
   void OnConBounds(
      int                index,              ///< AMPL index of constraint
      double             lb,                 ///< constraint left-hand-side
      double             ub                  ///< constraint right-hand-side
      )
   {
      assert(index >= 0);
      assert(index < probdata->nconss);

      // nonlinear constraints are first
      if( index < (int)nlconslin.size() )
      {
         if( !SCIPisInfinity(scip, -lb) )
         {
            SCIP_CALL_THROW( SCIPchgLhsNonlinear(scip, probdata->conss[index], lb) );
         }
         if( !SCIPisInfinity(scip,  ub) )
         {
            SCIP_CALL_THROW( SCIPchgRhsNonlinear(scip, probdata->conss[index], ub) );
         }
      }
      else
      {
         if( !SCIPisInfinity(scip, -lb) )
         {
            SCIP_CALL_THROW( SCIPchgLhsLinear(scip, probdata->conss[index], lb) );
         }
         if( !SCIPisInfinity(scip,  ub) )
         {
            SCIP_CALL_THROW( SCIPchgRhsLinear(scip, probdata->conss[index], ub) );
         }
      }
   }

   /// receive notification of the initial value for a variable
   void OnInitialValue(
      int                var_index,          ///< AMPL index of variable
      double             value               ///< initial primal value of variable
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
      int                /* con_index */,    ///< AMPL index of constraint
      double             /* value */         ///< initial dual value of constraint
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

   /// handling of suffices for variable and constraint flags and SOS constraints
   ///
   /// regarding SOS in AMPL, see https://ampl.com/faqs/how-can-i-use-the-solvers-special-ordered-sets-feature/
   /// we pass the .ref suffix as weight to the SOS constraint handlers
   /// for a SOS2, the weights determine the order of variables in the set
   template<typename T> class SuffixHandler
   {
   private:
      AMPLProblemHandler& amplph;

      // type of suffix that is handled, or IGNORE if unsupported suffix
      enum
      {
         IGNORE,
         CONSINITIAL,
         CONSSEPARATE,
         CONSENFORCE,
         CONSCHECK,
         CONSPROPAGATE,
         CONSDYNAMIC,
         CONSREMOVABLE,
         VARINITIAL,
         VARREMOVABLE,
         VARSOSNO,
         VARREF,
      } suffix;

   public:
      /// constructor
      SuffixHandler(
         AMPLProblemHandler& amplph_,        ///< problem handler
         fmt::StringRef      name,           ///< name of suffix
         mp::suf::Kind       kind            ///< whether suffix applies to var, cons, etc
         )
      : amplph(amplph_),
        suffix(IGNORE)
      {
         switch( kind )
         {
            case mp::suf::Kind::CON:
               if( strncmp(name.data(), "initial", name.size()) == 0 )
               {
                  suffix = CONSINITIAL;
               }
               else if( strncmp(name.data(), "separate", name.size()) == 0 )
               {
                  suffix = CONSSEPARATE;
               }
               else if( strncmp(name.data(), "enforce", name.size()) == 0 )
               {
                  suffix = CONSENFORCE;
               }
               else if( strncmp(name.data(), "check", name.size()) == 0 )
               {
                  suffix = CONSCHECK;
               }
               else if( strncmp(name.data(), "propagate", name.size()) == 0 )
               {
                  suffix = CONSPROPAGATE;
               }
               else if( strncmp(name.data(), "dynamic", name.size()) == 0 )
               {
                  suffix = CONSDYNAMIC;
               }
               else if( strncmp(name.data(), "removable", name.size()) == 0 )
               {
                  suffix = CONSREMOVABLE;
               }
               else
               {
                  SCIPverbMessage(amplph.scip, SCIP_VERBLEVEL_HIGH, NULL, "Unknown constraint suffix <%.*s>. Ignoring.\n", (int)name.size(), name.data());
               }
               break;

            case mp::suf::Kind::CON_BIT:
               SCIPverbMessage(amplph.scip, SCIP_VERBLEVEL_HIGH, NULL, "Unknown constraint bit suffix <%.*s>. Ignoring.\n", (int)name.size(), name.data());
               break;

            case mp::suf::Kind::VAR:
            {
               if( strncmp(name.data(), "initial", name.size()) == 0 )
               {
                  suffix = VARINITIAL;
               }
               else if( strncmp(name.data(), "removable", name.size()) == 0 )
               {
                  suffix = VARREMOVABLE;
               }
               else if( strncmp(name.data(), "sosno", name.size()) == 0 )
               {
                  // SOS membership
                  suffix = VARSOSNO;
               }
               else if( strncmp(name.data(), "ref", name.size()) == 0 )
               {
                  // SOS weights
                  suffix = VARREF;
                  amplph.sosweights.resize(amplph.probdata->nvars, 0);
               }
               else
               {
                  SCIPverbMessage(amplph.scip, SCIP_VERBLEVEL_HIGH, NULL, "Unknown variable suffix <%.*s>. Ignoring.\n", (int)name.size(), name.data());
               }
               break;

            case mp::suf::Kind::VAR_BIT:
               SCIPverbMessage(amplph.scip, SCIP_VERBLEVEL_HIGH, NULL, "Unknown variable bit suffix <%.*s>. Ignoring.\n", (int)name.size(), name.data());
               break;

            case mp::suf::Kind::OBJ:
               SCIPverbMessage(amplph.scip, SCIP_VERBLEVEL_HIGH, NULL, "Unknown objective suffix <%.*s>. Ignoring.\n", (int)name.size(), name.data());
               break;

            case mp::suf::Kind::OBJ_BIT:
               SCIPverbMessage(amplph.scip, SCIP_VERBLEVEL_HIGH, NULL, "Unknown objective bit suffix <%.*s>. Ignoring.\n", (int)name.size(), name.data());
               break;

            case mp::suf::Kind::PROBLEM:
               SCIPverbMessage(amplph.scip, SCIP_VERBLEVEL_HIGH, NULL, "Unknown problem suffix <%.*s>. Ignoring.\n", (int)name.size(), name.data());
               break;

            case mp::suf::Kind::PROB_BIT:
               SCIPverbMessage(amplph.scip, SCIP_VERBLEVEL_HIGH, NULL, "Unknown problem bit suffix <%.*s>. Ignoring.\n", (int)name.size(), name.data());
               break;
            }
         }
      }

      void SetValue(
         int             index,              ///< index of variable, constraint, etc
         T               value               ///< value of suffix
      )
      {
         assert(index >= 0);
         switch( suffix )
         {
            case IGNORE :
               return;

            case CONSINITIAL:
               SCIP_CALL_THROW( SCIPsetConsInitial(amplph.scip, amplph.probdata->conss[index], value == 1) );
               break;

            case CONSSEPARATE:
               SCIP_CALL_THROW( SCIPsetConsSeparated(amplph.scip, amplph.probdata->conss[index], value == 1) );
               break;

            case CONSENFORCE:
               SCIP_CALL_THROW( SCIPsetConsEnforced(amplph.scip, amplph.probdata->conss[index], value == 1) );
               break;

            case CONSCHECK:
               SCIP_CALL_THROW( SCIPsetConsChecked(amplph.scip, amplph.probdata->conss[index], value == 1) );
               break;

            case CONSPROPAGATE:
               SCIP_CALL_THROW( SCIPsetConsPropagated(amplph.scip, amplph.probdata->conss[index], value == 1) );
               break;

            case CONSDYNAMIC:
               SCIP_CALL_THROW( SCIPsetConsDynamic(amplph.scip, amplph.probdata->conss[index], value == 1) );
               break;

            case CONSREMOVABLE:
               SCIP_CALL_THROW( SCIPsetConsRemovable(amplph.scip, amplph.probdata->conss[index], value == 1) );
               break;

            case VARINITIAL:
               assert(index < amplph.probdata->nvars);
               SCIP_CALL_THROW( SCIPvarSetInitial(amplph.probdata->vars[index], value == 1) );
               break;

            case VARREMOVABLE:
               assert(index < amplph.probdata->nvars);
               SCIP_CALL_THROW( SCIPvarSetRemovable(amplph.probdata->vars[index], value == 1) );
               break;

            case VARSOSNO:
               // remember that variable index belongs to SOS identified by value
               amplph.sosvars[(int)value].push_back(index);
               break;

            case VARREF:
               // remember that variable index has weight value
               amplph.sosweights[index] = (int)value;
               break;
         }
      }
   };

   typedef SuffixHandler<int> IntSuffixHandler;
   /// receive notification of an integer suffix
   IntSuffixHandler OnIntSuffix(
      fmt::StringRef     name,               ///< suffix name, not null-terminated
      mp::suf::Kind      kind,               ///< suffix kind
      int                /*num_values*/      ///< number of values to expect
      )
   {
      return IntSuffixHandler(*this, name, kind);
   }

   typedef SuffixHandler<SCIP_Real> DblSuffixHandler;
   /// receive notification of a double suffix
   DblSuffixHandler OnDblSuffix(
      fmt::StringRef     name,               ///< suffix name, not null-terminated
      mp::suf::Kind      kind,               ///< suffix kind
      int                /*num_values*/      ///< number of values to expect
      )
   {
      return DblSuffixHandler(*this, name, kind);
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
         AMPLProblemHandler& amplph_,        ///< problem handler
         int                 constraintIndex_///< constraint index
         )
      : amplph(amplph_),
        constraintIndex(constraintIndex_)
      {
         assert(constraintIndex_ >= 0);
         assert(constraintIndex_ < amplph.probdata->nconss);
      }

      // constructor for linear objective
      explicit LinearPartHandler(
         AMPLProblemHandler& amplph_         ///< problem handler
         )
      : amplph(amplph_),
        constraintIndex(-1)
      { }

      void AddTerm(
         int             variableIndex,      ///< AMPL index of variable
         double          coefficient         ///< coefficient of variable
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
            SCIP_CONS* lincons = amplph.probdata->conss[constraintIndex];
            SCIP_CALL_THROW( SCIPaddCoefLinear(amplph.scip, lincons, amplph.probdata->vars[variableIndex], coefficient) );
         }
      }
   };

   typedef LinearPartHandler LinearObjHandler;

   /// receive notification of the linear part of an objective
   LinearPartHandler OnLinearObjExpr(
      int                objectiveIndex,     ///< index of objective
      int                /* numLinearTerms *////< number of terms to expect
      )
   {
      if( objectiveIndex >= 1 )
         OnUnhandled("multiple objective functions");

      return LinearObjHandler(*this);
   }

   typedef LinearPartHandler LinearConHandler;

   /// receive notification of the linear part of a constraint
   LinearConHandler OnLinearConExpr(
      int                constraintIndex,    ///< index of constraint
      int                /* numLinearTerms *////< number of terms to expect
      )
   {
      return LinearConHandler(*this, constraintIndex);
   }

   /// receives notification of a `Boolean value <mp::expr::BOOL>`
   LogicalExpr OnBool(
      bool               value
      )
   {
      SCIP_EXPR* expr;

      SCIP_CALL_THROW( SCIPcreateExprValue(scip, &expr, value ? 1.0 : 0.0, NULL, NULL) );

      // remember that we have to release this expr
      exprstorelease.push_back(expr);

      return expr;
   }

   /// receives notification of a `logical not <mp::expr::NOT>`
   LogicalExpr OnNot(
      LogicalExpr        arg
      )
   {
      SCIP_EXPR* expr;
      SCIP_VAR* var;
      SCIP_Bool val;

      LogicalExprToVarVal(arg, var, val);
      if( var != NULL )
      {
         SCIP_CALL_THROW( SCIPgetNegatedVar(scip, var, &var) );
         SCIP_CALL_THROW( SCIPcreateExprVar(scip, &expr, var, NULL, NULL) );
      }
      else
      {
         SCIP_CALL_THROW( SCIPcreateExprValue(scip, &expr, val ? 1.0 : 0.0, NULL, NULL) );
      }

      // remember that we have to release this expr
      exprstorelease.push_back(expr);

      return expr;
   }

   /// receives notification of a `binary logical expression <mp::expr::FIRST_BINARY_LOGICAL>`
   LogicalExpr OnBinaryLogical(
      mp::expr::Kind     kind,
      LogicalExpr        lhs,
      LogicalExpr        rhs
      )
   {
      SCIP_VAR* lhsvar = NULL;
      SCIP_VAR* rhsvar = NULL;
      SCIP_Bool lhsval;
      SCIP_Bool rhsval;
      SCIP_EXPR* expr;

      assert(lhs != NULL);
      assert(rhs != NULL);

      LogicalExprToVarVal(lhs, lhsvar, lhsval);
      LogicalExprToVarVal(rhs, rhsvar, rhsval);

      switch( kind )
      {
         case mp::expr::OR:
         {
            if( lhsvar == NULL && rhsvar == NULL )
            {
               SCIP_CALL_THROW( SCIPcreateExprValue(scip, &expr, lhsval != 0.0 || rhsval != 0.0 ? 1.0 : 0.0, NULL, NULL) );
               exprstorelease.push_back(expr);
               break;
            }

            if( (lhsvar == NULL && lhsval != 0.0) || (rhsvar == NULL && rhsval != 0.0) )
            {
               /* nonzero or rhs == 1, lhs or nonzero == 1 */
               SCIP_CALL_THROW( SCIPcreateExprValue(scip, &expr, 1.0, NULL, NULL) );
               exprstorelease.push_back(expr);
               break;
            }

            if( lhsvar == NULL )
            {
               /* zero or rhs == rhs */
               assert(lhsval == 0.0);
               expr = rhs;
               break;
            }

            if( rhsvar == NULL )
            {
               /* lhs or zero == lhs */
               assert(rhsval == 0.0);
               expr = lhs;
               break;
            }

            /* create new resvar and constraint resvar = lhsvar or rhsvar */
            SCIP_VAR* vars[2];
            SCIP_VAR* resvar;
            SCIP_CONS* cons;

            std::string name = std::string("_logic") + std::to_string((long long)logiccount++);
            SCIP_CALL_THROW( SCIPcreateVarBasic(scip, &resvar, name.c_str(), 0.0, 1.0, 0.0, SCIP_VARTYPE_BINARY) );
            SCIP_CALL_THROW( SCIPaddVar(scip, resvar) );
            SCIP_CALL_THROW( SCIPcreateExprVar(scip, &expr, resvar, NULL, NULL) );
            exprstorelease.push_back(expr);

            vars[0] = lhsvar;
            vars[1] = rhsvar;
            name += "def";
            SCIP_CALL_THROW( SCIPcreateConsBasicOr(scip, &cons, name.c_str(), resvar, 2, vars) );
            SCIP_CALL_THROW( SCIPaddCons(scip, cons) );

            SCIP_CALL_THROW( SCIPreleaseVar(scip, &resvar) );
            SCIP_CALL_THROW( SCIPreleaseCons(scip, &cons) );

            break;
         }

         case mp::expr::AND:
         {
            if( lhsvar == NULL && rhsvar == NULL )
            {
               SCIP_CALL_THROW( SCIPcreateExprValue(scip, &expr, lhsval != 0.0 && rhsval != 0.0 ? 1.0 : 0.0, NULL, NULL) );
               exprstorelease.push_back(expr);
               break;
            }

            if( (lhsvar == NULL && lhsval == 0.0) || (rhsvar == NULL && rhsval == 0.0) )
            {
               /* zero and rhs == 0, lhs and zero == 0 */
               SCIP_CALL_THROW( SCIPcreateExprValue(scip, &expr, 0.0, NULL, NULL) );
               exprstorelease.push_back(expr);
               break;
            }

            if( lhsvar == NULL )
            {
               /* nonzero and rhs == rhs */
               assert(lhsval != 0.0);
               expr = rhs;
               break;
            }

            if( rhsvar == NULL )
            {
               /* lhs and nonzero == lhs */
               assert(rhsval != 0.0);
               expr = lhs;
               break;
            }

            /* create new resvar and constraint resvar = lhsvar and rhsvar */
            SCIP_VAR* vars[2];
            SCIP_VAR* resvar;
            SCIP_CONS* cons;

            std::string name = std::string("_logic") + std::to_string((long long)logiccount++);
            SCIP_CALL_THROW( SCIPcreateVarBasic(scip, &resvar, name.c_str(), 0.0, 1.0, 0.0, SCIP_VARTYPE_BINARY) );
            SCIP_CALL_THROW( SCIPaddVar(scip, resvar) );
            SCIP_CALL_THROW( SCIPcreateExprVar(scip, &expr, resvar, NULL, NULL) );
            exprstorelease.push_back(expr);

            vars[0] = lhsvar;
            vars[1] = rhsvar;
            name += "def";
            SCIP_CALL_THROW( SCIPcreateConsBasicAnd(scip, &cons, name.c_str(), resvar, 2, vars) );
            SCIP_CALL_THROW( SCIPaddCons(scip, cons) );

            SCIP_CALL_THROW( SCIPreleaseVar(scip, &resvar) );
            SCIP_CALL_THROW( SCIPreleaseCons(scip, &cons) );

            break;
         }

         case mp::expr::IFF:
         {
            // the IFF operator returns 1 if both operands are nonzero or both are zero and returns zero otherwise
            // so this is lhs == rhs
            if( lhsvar == NULL && rhsvar == NULL )
            {
               SCIP_CALL_THROW( SCIPcreateExprValue(scip, &expr, lhsval == rhsval ? 1.0 : 0.0, NULL, NULL) );
               exprstorelease.push_back(expr);
               break;
            }

            if( lhsvar == NULL )
            {
               std::swap(lhs, rhs);
               std::swap(lhsval, rhsval);
               std::swap(lhsvar, rhsvar);
            }
            assert(lhsvar != NULL);

            if( rhsvar == NULL )
            {
               // expression is lhsvar == true
               // so we return lhsvar or ~lhsvar
               if( rhsval == TRUE )
               {
                  expr = lhs;
               }
               else
               {
                  SCIP_CALL_THROW( SCIPgetNegatedVar(scip, lhsvar, &lhsvar) );
                  SCIP_CALL_THROW( SCIPcreateExprVar(scip, &expr, lhsvar, NULL, NULL) );
                  exprstorelease.push_back(expr);
               }
               break;
            }

            // expressions is lhsvar == rhsvar
            // we create a new variable auxvar and add a constraint xor(auxvar, lhsvar, rhsvar, TRUE)
            // to ensure auxvar = (lhsvar == rhsvar)
            SCIP_VAR* vars[3];
            SCIP_CONS* cons;
            std::string name = std::string("_logic") + std::to_string((long long)logiccount++);
            SCIP_CALL_THROW( SCIPcreateVarBasic(scip, &vars[0], name.c_str(), 0.0, 1.0, 0.0, SCIP_VARTYPE_BINARY) );
            SCIP_CALL_THROW( SCIPaddVar(scip, vars[0]) );
            SCIP_CALL_THROW( SCIPcreateExprVar(scip, &expr, vars[0], NULL, NULL) );
            exprstorelease.push_back(expr);

            vars[1] = lhsvar;
            vars[2] = rhsvar;
            name += "def";
            SCIP_CALL_THROW( SCIPcreateConsBasicXor(scip, &cons, name.c_str(), TRUE, 3, vars) );
            SCIP_CALL_THROW( SCIPaddCons(scip, cons) );

            SCIP_CALL_THROW( SCIPreleaseVar(scip, &vars[0]) );
            SCIP_CALL_THROW( SCIPreleaseCons(scip, &cons) );

            break;
         }

         default:
            OnUnhandled(mp::expr::str(kind));
            return NULL;
      }

      return expr;
   }

   /// receives notification of a `relational expression <mp::expr::FIRST_RELATIONAL>`
   /// we only handle equality or inequality between binary variables and boolean values here
   LogicalExpr OnRelational(
      mp::expr::Kind     kind,
      NumericExpr        lhs,
      NumericExpr        rhs
      )
   {
      SCIP_VAR* lhsvar = NULL;
      SCIP_VAR* rhsvar = NULL;
      SCIP_Bool lhsval;
      SCIP_Bool rhsval;
      SCIP_EXPR* expr;

      assert(lhs != NULL);
      assert(rhs != NULL);

      LogicalExprToVarVal(lhs, lhsvar, lhsval);
      LogicalExprToVarVal(rhs, rhsvar, rhsval);

      switch( kind )
      {
         case mp::expr::EQ:
         case mp::expr::NE:
         {
            bool isne = (kind == mp::expr::NE);
            if( lhsvar == NULL && rhsvar == NULL )
            {
               SCIP_CALL_THROW( SCIPcreateExprValue(scip, &expr, lhsval == rhsval ? (isne ? 0.0 : 1.0) : (isne ? 1.0 : 0.0), NULL, NULL) );
               exprstorelease.push_back(expr);
               break;
            }

            if( lhsvar == NULL )
            {
               std::swap(lhs, rhs);
               std::swap(lhsval, rhsval);
               std::swap(lhsvar, rhsvar);
            }
            assert(lhsvar != NULL);

            if( rhsvar == NULL )
            {
               // expression is lhsvar == true or lhsvar == false if EQ
               // so we return lhsvar or ~lhsvar, opposite if NE
               if( rhsval == (isne ? FALSE : TRUE) )
               {
                  expr = lhs;
               }
               else
               {
                  SCIP_CALL_THROW( SCIPgetNegatedVar(scip, lhsvar, &lhsvar) );
                  SCIP_CALL_THROW( SCIPcreateExprVar(scip, &expr, lhsvar, NULL, NULL) );
                  exprstorelease.push_back(expr);
               }
               break;
            }

            // expressions is lhsvar == rhsvar or lhsvar != rhsvar
            // we create a new variable auxvar and add a constraint xor(auxvar, lhsvar, rhsvar, isne ? FALSE : TRUE)
            // to ensure auxvar = (lhsvar == rhsvar)  or  auxvar = (lhsvar != rhsvar)

            SCIP_VAR* vars[3];
            SCIP_CONS* cons;
            std::string name = std::string("_logic") + std::to_string((long long)logiccount++);
            SCIP_CALL_THROW( SCIPcreateVarBasic(scip, &vars[0], name.c_str(), 0.0, 1.0, 0.0, SCIP_VARTYPE_BINARY) );
            SCIP_CALL_THROW( SCIPaddVar(scip, vars[0]) );
            SCIP_CALL_THROW( SCIPcreateExprVar(scip, &expr, vars[0], NULL, NULL) );
            exprstorelease.push_back(expr);

            vars[1] = lhsvar;
            vars[2] = rhsvar;
            name += "def";
            SCIP_CALL_THROW( SCIPcreateConsBasicXor(scip, &cons, name.c_str(), isne ? FALSE : TRUE, 3, vars) );
            SCIP_CALL_THROW( SCIPaddCons(scip, cons) );

            SCIP_CALL_THROW( SCIPreleaseVar(scip, &vars[0]) );
            SCIP_CALL_THROW( SCIPreleaseCons(scip, &cons) );

            break;
         }

         default:
            OnUnhandled(mp::expr::str(kind));
            return NULL;
      }

      return expr;
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

         SCIP_CALL_THROW( SCIPcreateVarBasic(scip, &objvar, "nlobjvar", -SCIPinfinity(scip), SCIPinfinity(scip), 1.0, SCIP_VARTYPE_CONTINUOUS) );
         SCIP_CALL_THROW( SCIPaddVar(scip, objvar) );

         SCIP_CALL_THROW( SCIPcreateConsBasicNonlinear(scip, &objcons, "objcons", objexpr,
            SCIPgetObjsense(scip) == SCIP_OBJSENSE_MINIMIZE ? -SCIPinfinity(scip) : 0.0,
            SCIPgetObjsense(scip) == SCIP_OBJSENSE_MAXIMIZE ?  SCIPinfinity(scip) : 0.0) );
         SCIP_CALL_THROW( SCIPaddLinearVarNonlinear(scip, objcons, objvar, -1.0) );
         SCIP_CALL_THROW( SCIPaddCons(scip, objcons) );

         if( initsol != NULL )
         {
            /* compute value for objvar in initial solution from other variable values */
            SCIP_CALL_THROW( SCIPevalExpr(scip, objexpr, initsol, 0) );
            if( SCIPexprGetEvalValue(objexpr) != SCIP_INVALID )
            {
               SCIPsetSolVal(scip, initsol, objvar, SCIPexprGetEvalValue(objexpr));
            }
            else
            {
               SCIPwarningMessage(scip, "Objective function could not be evaluated in initial point. Domain error.");
            }
         }

         SCIP_CALL_THROW( SCIPreleaseCons(scip, &objcons) );
         SCIP_CALL_THROW( SCIPreleaseVar(scip, &objvar) );
      }

      // add linear terms to expressions of nonlinear constraints (should be ok to do this one-by-one for now)
      for( size_t i = 0; i < nlconslin.size(); ++i )
      {
         for( size_t j = 0; j < nlconslin[i].size(); ++j )
         {
            SCIP_CALL_THROW( SCIPaddLinearVarNonlinear(scip, probdata->conss[i], nlconslin[i][j].second, nlconslin[i][j].first) );
         }
      }

      // add constraints
      for( int i = 0; i < probdata->nconss; ++i )
      {
         SCIP_CALL_THROW( SCIPaddCons(scip, probdata->conss[i]) );
      }

      // add SOS constraints
      std::vector<SCIP_VAR*> setvars;     // variables in one SOS
      std::vector<SCIP_Real> setweights;  // weights for one SOS
      if( !sosvars.empty() )
      {
         setvars.resize(probdata->nvars);
         probdata->islp = false;
      }
      if( !sosweights.empty() )
         setweights.resize(probdata->nvars);
      for( std::map<int, std::vector<int> >::iterator sosit(sosvars.begin()); sosit != sosvars.end(); ++sosit )
      {
         assert(sosit->first != 0);
         assert(!sosit->second.empty());

         // a negative SOS identifier means SOS2
         bool issos2 = sosit->first < 0;

         if( issos2 && sosweights.empty() )
         {
            // if no .ref suffix was given for a SOS2 constraint, then we consider this as an error
            // since the weights determine the order
            // for a SOS1, the weights only specify branching preference, so can treat them as optional
            OnUnhandled("SOS2 requires variable .ref suffix");
         }

         for( size_t i = 0; i < sosit->second.size(); ++i )
         {
            int varidx = sosit->second[i];
            setvars[i] = probdata->vars[varidx];

            if( issos2 && sosweights[varidx] == 0 )
               // 0 is the default if no ref was given for a variable; we don't allow this for SOS2
               OnUnhandled("Missing .ref value for SOS2 variable");
            if( !sosweights.empty() )
               setweights[i] = (SCIP_Real)sosweights[varidx];
         }

         SCIP_CONS* cons;
         char name[20];
         if( !issos2 )
         {
            (void) SCIPsnprintf(name, 20, "sos1_%d", sosit->first);
            SCIP_CALL_THROW( SCIPcreateConsBasicSOS1(scip, &cons, name, sosit->second.size(), setvars.data(), setweights.empty() ? NULL : setweights.data()) );
         }
         else
         {
            (void) SCIPsnprintf(name, 20, "sos2_%d", -sosit->first);
            SCIP_CALL_THROW( SCIPcreateConsBasicSOS2(scip, &cons, name, sosit->second.size(), setvars.data(), setweights.data()) );
         }
         SCIP_CALL_THROW( SCIPaddCons(scip, cons) );
         SCIP_CALL_THROW( SCIPreleaseCons(scip, &cons) );
      }

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

class SCIPNLFeeder : public mp::NLFeeder<SCIPNLFeeder, SCIP_EXPR*>
{
private:
   SCIP*                 scip;               ///< SCIP data structure (problem to write)
   bool                  nlcomments;         ///< whether to write nl files with comments
   const char*           probname;           ///< problem name
   SCIP_Bool             transformed;        ///< is transformed problem written
   SCIP_OBJSENSE         objsense;           ///< objective sense
   SCIP_Real             objscale;           ///< objective scale
   SCIP_Real             objoffset;          ///< objective offset
   SCIP_VAR**            activevars;         ///< active variables
   int                   nactivevars;        ///< number of active variables
   SCIP_VAR**            fixedvars;          ///< fixed variables
   int                   nfixedvars;         ///< number of fixed variables
   SCIP_CONS**           allconss;           ///< all constraints given to writer
   int                   nallconss;          ///< number of all constraints
   SCIP_Bool             genericnames;       ///< are generic names used

   SCIP_CONSHDLR*        conshdlr_nonlinear; ///< nonlinear constraints handler
   SCIP_CONSHDLR*        conshdlr_linear;    ///< linear constraints handler
   SCIP_CONSHDLR*        conshdlr_setppc;    ///< setppc constraints handler
   SCIP_CONSHDLR*        conshdlr_logicor;   ///< logicor constraints handler
   SCIP_CONSHDLR*        conshdlr_knapsack;  ///< knapsack constraints handler
   SCIP_CONSHDLR*        conshdlr_varbound;  ///< varbound constraints handlers

   mp::NLHeader          nlheader;           ///< NL header with various counts
   SCIP_VAR**            vars;               ///< variables in AMPL order
   int                   nvars;              ///< number of variables (= nactivevars)
   SCIP_HASHMAP*         var2idx;            ///< map variable to AMPL index
   SCIP_CONS**           algconss;           ///< algebraic constraints that will be written, permuted in AMPL order
   SCIP_Real*            algconsslhs;        ///< left hand side of algebraic constraints
   SCIP_Real*            algconssrhs;        ///< right hand side of algebraic constraints
   int                   nalgconss;          ///< number of algebraic constraint we will actually write
   SCIP_VAR**            aggconss;           ///< fixed variable for which aggregation constraints need to be written
   int                   naggconss;          ///< number of fixed variables for which aggregation constraints are written

   /** variable types by which variables need to be ordered for .nl
    * (names are taken from pyomo nl writer, with those for nonlinear objective removed)
    */
   typedef enum
   {
      ConNonlinearVars = 0,     /* only in cons */
      ConNonlinearVarsInt = 1,  /* only in cons */
      LinearVars = 2,
      LinearVarsBool = 3,
      LinearVarsInt = 4
   } NlVarType;

   /** checks variable types and other properties for nlheader;
    * sets up variables permutation
    */
   void analyseVariables()
   {
      NlVarType* vartype = NULL;
      SCIP_HASHMAP* var2expr = NULL;

      int nlvars_cons = 0;
      int binvars_lin = 0;
      int intvars_lin = 0;
      int discrvars_nlcons = 0;

      nlheader.max_var_name_len = 0;

      /* number of nonzeros in objective gradient */
      nlheader.num_obj_nonzeros = 0;

      if( conshdlr_nonlinear != NULL )
         var2expr = SCIPgetVarExprHashmapNonlinear(conshdlr_nonlinear);

      SCIP_CALL_THROW( SCIPallocBufferArray(scip, &vartype, nactivevars + nfixedvars) );

      /* collect statistics on variables; determine variable types */
      for( int i = 0; i < nactivevars + nfixedvars; ++i )
      {
         SCIP_VAR* var = (i < nactivevars ? activevars[i] : fixedvars[i-nactivevars]);
         SCIP_Bool isdiscrete;
         SCIP_Bool isnonlinear;

         if( SCIPvarGetObj(var) != 0.0 )
            ++nlheader.num_obj_nonzeros;

         isdiscrete = SCIPvarGetType(var) <= SCIP_VARTYPE_INTEGER;

         /* we think of a variable as nonlinear if cons_nonlinear has a SCIP_EXPR* for this variable
          * this is usually an overestimation, since also variables that appear only linearly in nonlinear constraints
          * are regarded as nonlinear this way
          */
         isnonlinear = var2expr != NULL && SCIPhashmapExists(var2expr, (void*)var);

         /* this is how Pyomo counts vars (nlvars_* = nlvb,c,o) when writing NL
          * https://github.com/Pyomo/pyomo/blob/main/pyomo/repn/plugins/ampl/ampl_.py#L1202
          * this, together with the ominous line below, seems to correspond to what AMPL writes
          */
         if( isnonlinear )
         {
            /* nonlinear (in constraints only, as this is SCIP) */
            ++nlvars_cons;
            if( isdiscrete )
            {
               ++discrvars_nlcons;
               vartype[i] = ConNonlinearVarsInt;
            }
            else
               vartype[i] = ConNonlinearVars;
         }
         else
         {
            /* linear */
            if( isdiscrete )
            {
               /* for compatibility with AMPL generated nl files, count integer with 0/1 bounds as binary, too */
               if( !SCIPisNegative(scip, SCIPvarGetLbGlobal(var)) && SCIPisLE(scip, SCIPvarGetUbGlobal(var), 1.0) )
               {
                  ++binvars_lin;
                  vartype[i] = LinearVarsBool;
               }
               else
               {
                  ++intvars_lin;
                  vartype[i] = LinearVarsInt;
               }
            }
            else
               vartype[i] = LinearVars;
         }

         if( !genericnames )
         {
            int namelen = (int)strlen(SCIPvarGetName(var));
            if( namelen > nlheader.max_var_name_len )
               nlheader.max_var_name_len = namelen;
         }
      }

      /* setup var permutation */
      assert(vars == NULL);
      SCIP_CALL_THROW( SCIPallocBlockMemoryArray(scip, &vars, nactivevars + nfixedvars) );
      SCIP_CALL_THROW( SCIPhashmapCreate(&var2idx, SCIPblkmem(scip), nactivevars + nfixedvars) );
      nvars = 0;
      for( int vtype = ConNonlinearVars; vtype <= LinearVarsInt; ++vtype )
         for( int i = 0; i < nactivevars + nfixedvars; ++i )
            if( vartype[i] == (NlVarType)vtype )
            {
               vars[nvars] = (i < nactivevars ? activevars[i] : fixedvars[i-nactivevars]);
               SCIP_CALL_THROW( SCIPhashmapInsertInt(var2idx, (void*)vars[nvars], nvars) );
               ++nvars;
            }
      assert(nvars == nactivevars + nfixedvars);

      SCIPfreeBufferArray(scip, &vartype);

      nlheader.num_vars = nvars;

      /* number of nonlinear variables
       * setting num_nl_vars_in_objs = nlvars_cons looks odd, but makes the generated nl files
       * consistent with what AMPL or Pyomo writes
       */
      nlheader.num_nl_vars_in_cons = nlvars_cons;
      nlheader.num_nl_vars_in_objs = nlvars_cons;
      nlheader.num_nl_vars_in_both = 0;

      /* number of linear network variables */
      nlheader.num_linear_net_vars = 0;

      /* number of linear binary and integer variables */
      nlheader.num_linear_binary_vars = binvars_lin;
      nlheader.num_linear_integer_vars = intvars_lin;

      /* number of integer nonlinear variables */
      nlheader.num_nl_integer_vars_in_both = 0;
      nlheader.num_nl_integer_vars_in_cons = discrvars_nlcons;
      nlheader.num_nl_integer_vars_in_objs = 0;
   }

   /** checks constraint types and other properties for nlheader;
    * sets up constraints permutation
    */
   void analyzeConstraints()
   {
      /* collect algebraic constraints and their side: for AMPL, nonlinear comes before linear */
      SCIP_CALL_THROW( SCIPallocBlockMemoryArray(scip, &algconss, nallconss) );
      SCIP_CALL_THROW( SCIPallocBlockMemoryArray(scip, &algconsslhs, nallconss) );
      SCIP_CALL_THROW( SCIPallocBlockMemoryArray(scip, &algconssrhs, nallconss) );

      nalgconss = 0;
      if( nlheader.num_nl_vars_in_cons > 0 )
      {
         for( int i = 0; i < nallconss; ++i )
         {
            SCIP_CONS* cons = allconss[i];
            if( SCIPconsGetHdlr(cons) == conshdlr_nonlinear )
            {
               algconss[nalgconss] = cons;
               algconsslhs[nalgconss] = SCIPgetLhsNonlinear(cons);
               algconssrhs[nalgconss] = SCIPgetRhsNonlinear(cons);
               ++nalgconss;
            }
         }
      }
      /* total number of nonlinear constraints */
      nlheader.num_nl_cons = nalgconss;

      /* pick constraints we recognize as linear
       * count ranged and equality constraints
       * check constraint name lengths (if not skipped due to being generic)
       * count number of variables in constraints
       */
      nlheader.num_ranges = 0; /* number of ranged constraints */
      nlheader.num_eqns = 0;   /* number of equality constraints */
      nlheader.max_con_name_len = 0;  /* maximal length of constraints' names */
      nlheader.num_con_nonzeros = 0;  /* number of nonzeros in constraints' Jacobian */
      for( int i = 0; i < nallconss; ++i )
      {
         SCIP_CONS* cons = allconss[i];
         SCIP_CONSHDLR* conshdlr = SCIPconsGetHdlr(cons);
         SCIP_Real lhs;
         SCIP_Real rhs;

         if( conshdlr == conshdlr_nonlinear )
         {
            lhs = SCIPgetLhsNonlinear(cons);
            rhs = SCIPgetRhsNonlinear(cons);
         }
         else
         {
            /* negated variables in may not show up in fixedvars
             * so we instead replace the negation when providing the coefficients of the linear constraint
             * this means additional constants to subtract from lhs/rhs
             */
            if( conshdlr == conshdlr_linear )
            {
               int nconsvars = SCIPgetNVarsLinear(scip, cons);
               SCIP_VAR** consvars = SCIPgetVarsLinear(scip, cons);
               SCIP_Real* conscoefs = SCIPgetValsLinear(scip, cons);
               SCIP_Real negconstant = 0.0;
               for( int v = 0; v < nconsvars; ++v )
                  if( SCIPvarIsNegated(consvars[v]) )
                     negconstant += conscoefs[v];

               lhs = SCIPgetLhsLinear(scip, cons);
               if( !SCIPisInfinity(scip, -lhs) )
                  lhs -= negconstant;

               rhs = SCIPgetRhsLinear(scip, cons);
               if( !SCIPisInfinity(scip, rhs) )
                  rhs -= negconstant;
            }
            else if( conshdlr == conshdlr_setppc )
            {
               int nconsvars = SCIPgetNVarsSetppc(scip, cons);
               SCIP_VAR** consvars = SCIPgetVarsSetppc(scip, cons);
               SCIP_Real negconstant = 0.0;
               for( int v = 0; v < nconsvars; ++v )
                  if( SCIPvarIsNegated(consvars[v]) )
                     negconstant += 1.0;

               switch( SCIPgetTypeSetppc(scip, cons) )
               {
                  case SCIP_SETPPCTYPE_PARTITIONING:
                     lhs = 1.0 - negconstant;
                     rhs = 1.0 - negconstant;
                     break;
                  case SCIP_SETPPCTYPE_COVERING:
                     lhs = 1.0 - negconstant;
                     rhs = SCIPinfinity(scip);
                     break;
                  case SCIP_SETPPCTYPE_PACKING:
                     lhs = -SCIPinfinity(scip);
                     rhs = 1.0 - negconstant;
                     break;
               }
            }
            else if( conshdlr == conshdlr_logicor )
            {
               int nconsvars = SCIPgetNVarsLogicor(scip, cons);
               SCIP_VAR** consvars = SCIPgetVarsLogicor(scip, cons);
               SCIP_Real negconstant = 0.0;
               for( int v = 0; v < nconsvars; ++v )
                  if( SCIPvarIsNegated(consvars[v]) )
                     negconstant += 1.0;

               lhs = 1.0 - negconstant;
               rhs = SCIPinfinity(scip);
            }
            else if( conshdlr == conshdlr_knapsack )
            {
               int nconsvars = SCIPgetNVarsKnapsack(scip, cons);
               SCIP_VAR** consvars = SCIPgetVarsKnapsack(scip, cons);
               SCIP_Longint* weights = SCIPgetWeightsKnapsack(scip, cons);
               SCIP_Longint negweights = 0.0;
               for( int v = 0; v < nconsvars; ++v )
                  if( SCIPvarIsNegated(consvars[v]) )
                     negweights += weights[v];

               lhs = -SCIPinfinity(scip);
               rhs = (SCIP_Real)(SCIPgetCapacityKnapsack(scip, cons) - negweights);
            }
            else if( conshdlr == conshdlr_varbound )
            {
               /* lhs <= var + vbdcoef*vbdvar <= rhs */
               SCIP_Real negconstant = 0.0;
               if( SCIPvarIsNegated(SCIPgetVarVarbound(scip, cons)) )
                  negconstant = 1.0;
               if( SCIPvarIsNegated(SCIPgetVbdvarVarbound(scip, cons)) )
                  negconstant += SCIPgetVbdcoefVarbound(scip, cons);

               lhs = SCIPgetLhsVarbound(scip, cons);
               if( !SCIPisInfinity(scip, -lhs) )
                  lhs -= negconstant;

               rhs = SCIPgetRhsVarbound(scip, cons);
               if( !SCIPisInfinity(scip, rhs) )
                  rhs -= negconstant;
            }
            else
            {
               SCIPwarningMessage(scip, "constraint <%s> of type <%s> cannot be printed in requested format\n", SCIPconsGetName(cons), SCIPconshdlrGetName(conshdlr));
               continue;
            }
            algconss[nalgconss] = cons;
            algconsslhs[nalgconss] = lhs;
            algconssrhs[nalgconss] = rhs;
            ++nalgconss;
         }

         if( !SCIPisInfinity(scip, -lhs) && !SCIPisInfinity(scip, rhs) )
         {
            if( SCIPisEQ(scip, lhs, rhs) )
               ++nlheader.num_eqns;
            else
               ++nlheader.num_ranges;
         }

         if( !genericnames )
         {
            int namelen = (int)strlen(SCIPconsGetName(allconss[i]));
            if( namelen > nlheader.max_con_name_len )
               nlheader.max_con_name_len = namelen;
         }

         SCIP_Bool success;
         int nvarsincons;
         SCIP_CALL_THROW( SCIPgetConsNVars(scip, cons, &nvarsincons, &success) );
         if( !success )
         {
            /* this should never happen */
            SCIPwarningMessage(scip, "could not get number of variable from constraint handler <%s>; nonzero count in nl file will be wrong\n", SCIPconshdlrGetName(conshdlr));
         }
         else
         {
            nlheader.num_con_nonzeros += nvarsincons;
         }
      }
      assert(nalgconss <= nallconss);

      /* now add counts for aggregation constraints (definition of fixedvars that are aggregated, multiaggregated, or negated) */
      SCIP_CALL_THROW( SCIPallocBlockMemoryArray(scip, &aggconss, nfixedvars) );
      naggconss = 0;
      for( int i = 0; i < nfixedvars; ++i )
      {
         SCIP_VAR* var = fixedvars[i];

         switch( SCIPvarGetStatus(var) )
         {
            case SCIP_VARSTATUS_FIXED:
               continue;

            case SCIP_VARSTATUS_AGGREGATED:
            case SCIP_VARSTATUS_NEGATED:
               nlheader.num_con_nonzeros += 2;
               break;

            case SCIP_VARSTATUS_MULTAGGR:
               nlheader.num_con_nonzeros += SCIPvarGetMultaggrNVars(var) + 1;
               break;

            default:
               SCIPerrorMessage("unexpected variable status %d of fixed variable <%s>\n", SCIPvarGetStatus(var), SCIPvarGetName(var));
               SCIP_CALL_THROW( SCIP_ERROR );
         }

         if( !genericnames )
         {
            // AMPL constraint will be named aggr_<varname>
            int namelen = (int)strlen(SCIPvarGetName(var)) + 5;
            if( namelen > nlheader.max_con_name_len )
               nlheader.max_con_name_len = namelen;
         }

         aggconss[naggconss] = var;
         ++naggconss;

         ++nlheader.num_eqns;
      }

      nlheader.num_algebraic_cons = nalgconss + naggconss;
      nlheader.num_logical_cons = 0;

      /* no complementarity conditions */
      nlheader.num_compl_conds = 0;
      nlheader.num_nl_compl_conds = 0;
      nlheader.num_compl_dbl_ineqs = 0;
      nlheader.num_compl_vars_with_nz_lb = 0;

      /** no network constraints */
      nlheader.num_nl_net_cons = 0;
      nlheader.num_linear_net_cons = 0;
   }

   /* gets AMPL index of variable (using var2idx) */
   int getVarAMPLIndex(
      SCIP_VAR*          var
      )
   {
      int varidx = SCIPhashmapGetImageInt(var2idx, (void*)var);
      assert(varidx >= 0);
      assert(varidx != INT_MAX);
      assert(varidx < nvars);
      assert(vars[varidx] == var);
      return varidx;
   }

public:
   /// Constructor
   SCIPNLFeeder(
      SCIP*              scip_,              ///< SCIP data structure
      SCIP_Bool          nlbinary_,          ///< whether to write binary or text nl
      SCIP_Bool          nlcomments_,        ///< whether to include comments into nl
      const char*        probname_,          ///< problem name
      SCIP_Bool          transformed_,       ///< is transformed problem written
      SCIP_OBJSENSE      objsense_,          ///< objective sense
      SCIP_Real          objscale_,          ///< objective scale
      SCIP_Real          objoffset_,         ///< objective offset
      SCIP_VAR**         vars_,              ///< active variables
      int                nvars_,             ///< number of active variables
      SCIP_VAR**         fixedvars_,         ///< fixed variables
      int                nfixedvars_,        ///< number of fixed variables
      SCIP_CONS**        conss_,             ///< constraints
      int                nconss_,            ///< number of constraints
      SCIP_Bool          genericnames_       ///< are generic names used
   )
   : scip(scip_),
     nlcomments(nlcomments_),
     probname(probname_),
     transformed(transformed_),
     objsense(objsense_),
     objscale(objscale_),
     objoffset(objoffset_),
     activevars(vars_),
     nactivevars(nvars_),
     fixedvars(fixedvars_),
     nfixedvars(nfixedvars_),
     allconss(conss_),
     nallconss(nconss_),
     genericnames(genericnames_),
     vars(NULL),
     nvars(0),
     var2idx(NULL),
     algconss(NULL),
     algconsslhs(NULL),
     algconssrhs(NULL),
     nalgconss(0),
     aggconss(NULL),
     naggconss(0)
   {
      nlheader.format = nlbinary_ ? mp::NLHeader::BINARY : mp::NLHeader::TEXT;

      conshdlr_nonlinear = SCIPfindConshdlr(scip, "nonlinear");
      conshdlr_linear = SCIPfindConshdlr(scip, "linear");
      conshdlr_setppc = SCIPfindConshdlr(scip, "setppc");
      conshdlr_logicor = SCIPfindConshdlr(scip, "logicor");
      conshdlr_knapsack = SCIPfindConshdlr(scip, "knapsack");
      conshdlr_varbound = SCIPfindConshdlr(scip, "varbound");
   }

   ~SCIPNLFeeder()
   {
      SCIPfreeBlockMemoryArrayNull(scip, &aggconss, nfixedvars);
      SCIPfreeBlockMemoryArrayNull(scip, &algconssrhs, nallconss);
      SCIPfreeBlockMemoryArrayNull(scip, &algconsslhs, nallconss);
      SCIPfreeBlockMemoryArrayNull(scip, &algconss, nallconss);
      SCIPfreeBlockMemoryArrayNull(scip, &vars, nactivevars + nfixedvars);
      SCIPhashmapFree(&var2idx);
   }

   /** Provide NLHeader.
    *
    *  This method is called first.
    *
    *  NLHeader summarizes the model and provides some technical parameters,
    *  such as text/binary NL format.
    */
   mp::NLHeader Header()
   {
      analyseVariables();
      analyzeConstraints();

      nlheader.prob_name = probname;

      /* number of objectives
       * if objective is all zero in SCIP, then just don't write any objective to nl
       */
      if( nlheader.num_obj_nonzeros == 0 && objoffset == 0.0 )
         nlheader.num_objs = 0;
      else
         nlheader.num_objs = 1;
      nlheader.num_nl_objs = 0;

      /* number of functions */
      nlheader.num_funcs = 0;

      /* it would have been nice to handle fixed variables as common expressions,
       * but as common expression are handled like nonlinear expressions,
       * this would turn any linear constraint with fixed variables into common expressions
       */
      nlheader.num_common_exprs_in_both = 0;
      nlheader.num_common_exprs_in_cons = 0;
      nlheader.num_common_exprs_in_objs = 0;
      nlheader.num_common_exprs_in_single_cons = 0;
      nlheader.num_common_exprs_in_single_objs = 0;

      return nlheader;
   }

   /// NL comments?
   bool WantNLComments() const
   {
      return nlcomments;
   }

   /// currently we do not want to write size of each column in Jacobian
   /// (i.e., number of constraints each variable appears in)
   int WantColumnSizes() const
   {
      return 0;
   }

   int ObjType(
      int
   )
   {
      return SCIPgetObjsense(scip) == SCIP_OBJSENSE_MAXIMIZE ? 1 : 0;
   }

   template <class ObjGradWriter>
   void FeedObjGradient(
      int                i,
      ObjGradWriter&     gw
   )
   {
      assert(i == 0);
      auto gvw = gw.MakeVectorWriter(nlheader.num_obj_nonzeros);
      for( int v = 0; v < nvars; ++v )
      {
         SCIP_Real coef = SCIPvarGetObj(vars[v]);
         if( coef != 0.0 )
            gvw.Write(v, objscale * coef);
      }
   }

   template <class ObjExprWriter>
   void FeedObjExpression(
      int                i,
      ObjExprWriter&     ew
   )
   {
      assert(i == 0);
      ew.NPut(objscale * objoffset);
   }

   template <class VarBoundsWriter>
   void FeedVarBounds(
      VarBoundsWriter&   vbw
   )
   {
      for( int v = 0; v < nvars; ++v )
      {
         SCIP_Real lb = SCIPvarGetLbGlobal(vars[v]);
         SCIP_Real ub = SCIPvarGetUbGlobal(vars[v]);

         if( SCIPisInfinity(scip, -lb) )
            lb = -INFINITY;

         if( SCIPisInfinity(scip, ub) )
            ub = INFINITY;

         vbw.WriteLbUb(lb, ub);
      }
   }

   template <class ConBoundsWriter>
   void FeedConBounds(
      ConBoundsWriter& cbw
   )
   {
      for( int c = 0; c < nalgconss; ++c )
      {
         AlgConRange bnd;
         bnd.L = SCIPisInfinity(scip, -algconsslhs[c]) ? -INFINITY : algconsslhs[c];
         bnd.U = SCIPisInfinity(scip,  algconssrhs[c]) ?  INFINITY : algconssrhs[c];
         cbw.WriteAlgConRange(bnd);
      }

      for( int v = 0; v < naggconss; ++v )
      {
         SCIP_VAR* var = aggconss[v];
         AlgConRange bnd;

         switch( SCIPvarGetStatus(var) )
         {
            case SCIP_VARSTATUS_AGGREGATED:
               bnd.L = SCIPvarGetAggrConstant(var);
               break;

            case SCIP_VARSTATUS_NEGATED:
               bnd.L = SCIPvarGetNegationConstant(var);
               break;

            case SCIP_VARSTATUS_MULTAGGR:
               bnd.L = SCIPvarGetMultaggrConstant(var);
               break;

            default:
               SCIPerrorMessage("unexpected variable status %d of aggregated variable <%s>\n", SCIPvarGetStatus(var), SCIPvarGetName(var));
               SCIP_CALL_THROW( SCIP_ERROR );
         }

         bnd.U = bnd.L;
         cbw.WriteAlgConRange(bnd);
      }
   }

   /* this is for the comments in .nl files if comments enabled */
   const char* ConDescription(
      int i
   )
   {
      if( i < nalgconss )
         return SCIPconsGetName(algconss[i]);

      assert(i < nalgconss + naggconss);
      return SCIPvarGetName(aggconss[i-nalgconss]);
   }

   template <class ConLinearExprWriter>
   void FeedLinearConExpr(
      int                  i,
      ConLinearExprWriter& clw
   )
   {
      if( i < nlheader.num_nl_cons )
         return;

      if( i < nalgconss )
      {
         SCIP_CONS* cons = algconss[i];
         SCIP_CONSHDLR* conshdlr = SCIPconsGetHdlr(cons);

         if( conshdlr == conshdlr_linear )
         {
            SCIP_Real* conscoefs = SCIPgetValsLinear(scip, cons);
            SCIP_VAR** consvars = SCIPgetVarsLinear(scip, cons);
            int nconsvars = SCIPgetNVarsLinear(scip, cons);

            auto vw = clw.MakeVectorWriter(nconsvars);
            for( int v = 0; v < nconsvars; ++v )
               if( SCIPvarIsNegated(consvars[v]) )
                  vw.Write(getVarAMPLIndex(SCIPvarGetNegationVar(consvars[v])), -conscoefs[v]);
               else
                  vw.Write(getVarAMPLIndex(consvars[v]), conscoefs[v]);
         }
         else if( conshdlr == conshdlr_setppc )
         {
            SCIP_VAR** consvars = SCIPgetVarsSetppc(scip, cons);
            int nconsvars = SCIPgetNVarsSetppc(scip, cons);

            auto vw = clw.MakeVectorWriter(nconsvars);
            for( int v = 0; v < nconsvars; ++v )
               if( SCIPvarIsNegated(consvars[v]) )
                  vw.Write(getVarAMPLIndex(SCIPvarGetNegationVar(consvars[v])), -1.0);
               else
                  vw.Write(getVarAMPLIndex(consvars[v]), 1.0);
         }
         else if( conshdlr == conshdlr_logicor )
         {
            SCIP_VAR** consvars = SCIPgetVarsLogicor(scip, cons);
            int nconsvars = SCIPgetNVarsLogicor(scip, cons);

            auto vw = clw.MakeVectorWriter(nconsvars);
            for( int v = 0; v < nconsvars; ++v )
               if( SCIPvarIsNegated(consvars[v]) )
                  vw.Write(getVarAMPLIndex(SCIPvarGetNegationVar(consvars[v])), -1.0);
               else
                  vw.Write(getVarAMPLIndex(consvars[v]), 1.0);
         }
         else if( conshdlr == conshdlr_knapsack )
         {
            SCIP_Longint* weights = SCIPgetWeightsKnapsack(scip, cons);
            SCIP_VAR** consvars = SCIPgetVarsKnapsack(scip, cons);
            int nconsvars = SCIPgetNVarsKnapsack(scip, cons);

            auto vw = clw.MakeVectorWriter(nconsvars);
            for( int v = 0; v < nconsvars; ++v )
               if( SCIPvarIsNegated(consvars[v]) )
                  vw.Write(getVarAMPLIndex(SCIPvarGetNegationVar(consvars[v])), -(SCIP_Real)weights[v]);
               else
                  vw.Write(getVarAMPLIndex(consvars[v]), (SCIP_Real)weights[v]);
         }
         else
         {
            assert(conshdlr == conshdlr_varbound);

            auto vw = clw.MakeVectorWriter(2);
            if( SCIPvarIsNegated(SCIPgetVarVarbound(scip, cons)) )
               vw.Write(getVarAMPLIndex(SCIPvarGetNegationVar(SCIPgetVarVarbound(scip, cons))), -1.0);
            else
               vw.Write(getVarAMPLIndex(SCIPgetVarVarbound(scip, cons)), 1.0);

            if( SCIPvarIsNegated(SCIPgetVbdvarVarbound(scip, cons)) )
               vw.Write(getVarAMPLIndex(SCIPvarGetNegationVar(SCIPgetVbdvarVarbound(scip, cons))), -SCIPgetVbdcoefVarbound(scip, cons));
            else
               vw.Write(getVarAMPLIndex(SCIPgetVbdvarVarbound(scip, cons)), SCIPgetVbdcoefVarbound(scip, cons));
         }

         return;
      }

      assert(i < nalgconss + naggconss);
      SCIP_VAR* var = aggconss[i-nalgconss];

      switch( SCIPvarGetStatus(var) )
      {
         case SCIP_VARSTATUS_AGGREGATED:
         {
            /* var - aggrscalar*aggrvar = aggrconstant */
            auto vw = clw.MakeVectorWriter(2);
            vw.Write(getVarAMPLIndex(var), 1.0);
            vw.Write(getVarAMPLIndex(SCIPvarGetAggrVar(var)), -SCIPvarGetAggrScalar(var));
            break;
         }

         case SCIP_VARSTATUS_NEGATED:
         {
            /* var + negationvar = negationconstant */
            auto vw = clw.MakeVectorWriter(2);
            vw.Write(getVarAMPLIndex(var), 1.0);
            vw.Write(getVarAMPLIndex(SCIPvarGetNegationVar(var)), 1.0);
            break;
         }

         case SCIP_VARSTATUS_MULTAGGR:
         {
            /* var - sum_i aggrscalar_i aggrvar_i = aggrconstant */
            auto vw = clw.MakeVectorWriter(SCIPvarGetMultaggrNVars(var) + 1);
            vw.Write(getVarAMPLIndex(var), 1.0);
            for( int v = 0; v < SCIPvarGetMultaggrNVars(var); ++v )
               vw.Write(getVarAMPLIndex(SCIPvarGetMultaggrVars(var)[v]), -SCIPvarGetMultaggrScalars(var)[v]);
            break;
         }

         default:
         {
            SCIPerrorMessage("unexpected variable status %d of aggregated variable <%s>\n", SCIPvarGetStatus(var), SCIPvarGetName(var));
            SCIP_CALL_THROW( SCIP_ERROR );
         }
      }
   }

   template <class ConExprWriter>
   void FeedConExpression(
      int                i,
      ConExprWriter&     ew
   )
   {
      if( i >= nlheader.num_nl_cons )
      {
         ew.NPut(0.0);
         return;
      }

      SCIP_EXPR* rootexpr = SCIPgetExprNonlinear(algconss[i]);

      SCIP_EXPRITER* it;
      SCIP_CALL_THROW( SCIPcreateExpriter(scip, &it) );

      SCIP_CALL_THROW( SCIPexpriterInit(it, rootexpr, SCIP_EXPRITER_DFS, TRUE) );
      SCIPexpriterSetStagesDFS(it, SCIP_EXPRITER_ALLSTAGES);

      for( SCIP_EXPR* expr = SCIPexpriterGetCurrent(it); !SCIPexpriterIsEnd(it); expr = SCIPexpriterGetNext(it) )
      {
         switch( SCIPexpriterGetStageDFS(it) )
         {
            case SCIP_EXPRITER_ENTEREXPR:
            {
               // retrieve the ConExprWriter of parent expr
               ConExprWriter* parentew;
               if( expr == rootexpr )
                  parentew = &ew;
               else
                  parentew = (ConExprWriter*)SCIPexpriterGetExprUserData(it, SCIPexpriterGetParentDFS(it)).ptrval;
               assert(parentew != NULL);

               ConExprWriter* newew = NULL;

               if( SCIPisExprVar(scip, expr) )
               {
                  // TODO handle negated
                  SCIP_VAR* var = SCIPgetVarExprVar(expr);
                  parentew->VPut(getVarAMPLIndex(var), SCIPvarGetName(var));
               }
               else if( SCIPisExprValue(scip, expr) )
               {
                  parentew->NPut(SCIPgetValueExprValue(expr));
               }
               else if( SCIPisExprSum(scip, expr) )
               {
                  int nargs = SCIPexprGetNChildren(expr);

                  if( SCIPgetConstantExprSum(expr) != 0.0 )
                  {
                     if( nargs == 0 )
                     {
                        parentew->NPut(SCIPgetConstantExprSum(expr));
                        SCIP_EXPRITER_USERDATA userdata;
                        userdata.ptrval = NULL;
                        SCIPexpriterSetCurrentUserData(it, userdata);
                        break;
                     }

                     ++nargs;
                  }

                  // we will need to store two ConExprWriter's for sum or add
                  // one for the sum, and one for multiplication (coef*expr) of the currently considered child
                  // the one for the sum will go second, so in the child, we don't need to check for case of sum
                  // there is no default constructor for ConExprWriter, so we only alloc mem and then use replacement-new
                  SCIP_CALL_THROW( SCIPallocClearBufferArray(scip, &newew, 2) );

                  // this code is already complex enough for me that I skip some optimization for nargs==1 here
                  if( nargs == 2 )
                  {
                     new (newew+1) ConExprWriter(parentew->OPut2(mp::nl::ADD));
                  }
                  else
                  {
                     new (newew+1) ConExprWriter(parentew->OPutN(mp::nl::SUM, nargs));
                  }

                  if( SCIPgetConstantExprSum(expr) != 0.0 )
                     newew[1].NPut(SCIPgetConstantExprSum(expr));
               }
               else if( SCIPisExprProduct(scip, expr) )
               {
                  int nargs = SCIPexprGetNChildren(expr);

                  if( nargs == 2 )
                  {
                     newew = new ConExprWriter(parentew->OPut2(mp::nl::MUL));
                  }
                  else
                  {
                     //TODO
                     throw mp::UnsupportedError("can only do products of 2 factors at the moment");
                  }
               }
               else if( SCIPisExprPower(scip, expr) )
               {
                  if( SCIPgetExponentExprPow(expr) == 2.0 )
                     newew = new ConExprWriter(parentew->OPut1(mp::nl::POW2));
                  else
                     newew = new ConExprWriter(parentew->OPut2(mp::nl::POW_CONST_EXP));
               }
               else
               {
                  //TODO stop gracefully before throwing exception
                  throw mp::UnsupportedError("unsupported expr handler");
               }

               SCIP_EXPRITER_USERDATA userdata;
               userdata.ptrval = newew;
               SCIPexpriterSetCurrentUserData(it, userdata);

               break;
            }

            case SCIP_EXPRITER_VISITINGCHILD:
            {
               if( SCIPisExprSum(scip, expr) )
               {
                  int childidx = SCIPexpriterGetChildIdxDFS(it);
                  SCIP_Real coef = SCIPgetCoefsExprSum(expr)[childidx];

                  ConExprWriter* ews = (ConExprWriter*)SCIPexpriterGetCurrentUserData(it).ptrval;

                  if( coef != 1.0 )
                  {
                     // if coef, then create MUL and store ExprWriter in ews[0]
                     new (ews) ConExprWriter(ews[1].OPut2(mp::nl::MUL));
                     ews[0].NPut(coef);
                  }
                  else
                  {
                     // if trivial coef, then only move ews[1] (ExprWriter for SUM/ADD) to ews[0] (implementation forbids copy)
                     // cannot use move-assignment, because it asserts that destination and source have same nlw_, but my destination is not initialized
                     //ews[0] = std::move(ews[1]);
                     memcpy((void*)ews, (void*)(ews+1), sizeof(ConExprWriter));
                  }
               }
               break;
            }

            case SCIP_EXPRITER_VISITEDCHILD:
            {
               if( SCIPisExprSum(scip, expr) )
               {
                  int childidx = SCIPexpriterGetChildIdxDFS(it);

                  ConExprWriter* ews = (ConExprWriter*)SCIPexpriterGetCurrentUserData(it).ptrval;

                  if( SCIPgetCoefsExprSum(expr)[childidx] != 1.0 )
                  {
                     // destructor for ExprWriter that was stored for MUL
                     ews->~ConExprWriter();
                  }
                  else
                  {
                     // move ExprWrite for SUM back into 2nd position
                     //ews[1] = std::move(ews[0]);
                     memcpy((void*)(ews+1), (void*)ews, sizeof(ConExprWriter));
                  }
               }
               break;
            }

            case SCIP_EXPRITER_LEAVEEXPR:
            {
               ConExprWriter* ews = (ConExprWriter*)SCIPexpriterGetCurrentUserData(it).ptrval;
               if( SCIPisExprSum(scip, expr) )
               {
                  // destructor for ExprWriter for SUM/ADD
                  ews[1].~ConExprWriter();
                  SCIPfreeBufferArray(scip, &ews);
               }
               else
               {
                  // write exponent of power (if not 2)
                  if( SCIPisExprPower(scip, expr) && SCIPgetExponentExprPow(expr) != 2.0 )
                     ews->NPut(SCIPgetExponentExprPow(expr));

                  delete ews;
               }
               break;
            }
         }
      }

      SCIPfreeExpriter(&it);
   }

   template <class RowObjNameWriter>
   void FeedRowAndObjNames(
      RowObjNameWriter&  wrt
   )
   {
      if( !wrt || genericnames )
         return;

      for( int c = 0; c < nalgconss; ++c )
         wrt << SCIPconsGetName(algconss[c]);

      for( int v = 0; v < naggconss; ++v )
      {
         std::string aggname("aggr_");
         aggname += SCIPvarGetName(aggconss[v]);
         wrt << aggname.c_str();
      }

      wrt << "obj";
   }

   template <class ColNameWriter>
   void FeedColNames(
      ColNameWriter&     wrt
   )
   {
      if( !wrt || genericnames )
         return;

      for( int v = 0; v < nvars; ++v )
         wrt << SCIPvarGetName(vars[v]);
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
   assert((*probdata)->conss != NULL || (*probdata)->conss == 0);

   for( i = 0; i < (*probdata)->nconss; ++i )
   {
      SCIP_CALL( SCIPreleaseCons(scip, &(*probdata)->conss[i]) );
   }
   SCIPfreeBlockMemoryArrayNull(scip, &(*probdata)->conss, (*probdata)->nconss);

   for( i = 0; i < (*probdata)->nvars; ++i )
   {
      SCIP_CALL( SCIPreleaseVar(scip, &(*probdata)->vars[i]) );
   }
   SCIPfreeBlockMemoryArrayNull(scip, &(*probdata)->vars, (*probdata)->nvars);

   SCIPfreeBlockMemoryArrayNull(scip, &(*probdata)->filenamestub, (*probdata)->filenamestublen+5);

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

   *result = SCIP_DIDNOTRUN;

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
      catch( const mp::Error& e )
      {
         // some other error from ampl/mp, maybe invalid .nl file
         SCIPerrorMessage("%s\n", e.what());

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

/** problem writing method of reader */
static
SCIP_DECL_READERWRITE(readerWriteNl)
{  /*lint --e{715}*/
   mp::WriteNLResult writerresult;

   *result = SCIP_DIDNOTRUN;

   SCIPNLFeeder nlf(scip, FALSE, TRUE,
      name, transformed, objsense, objscale, objoffset,
      vars, nvars, fixedvars, nfixedvars,
      conss, nconss, genericnames);

   try
   {
      /* we need to give the NLWriter a filename, but got only a FILE* from SCIP
       * so we let the NLWriter write to a temporary file and then copy its content to file
       * TODO this way, the files with row and col names are lost
       */
      char tempname[L_tmpnam+4];
      FILE* tempfile;

      if( std::tmpnam(tempname) == NULL )
      {
         SCIPerrorMessage("Cannot generate name for temporary file: error %d\n", errno);
         return SCIP_FILECREATEERROR;
      }

      mp::NLUtils nlutils;
      writerresult = mp::WriteNLFile(tempname, nlf, nlutils);

      switch( writerresult.first )
      {
         case NLW2_WriteNL_OK:
            break;
         case NLW2_WriteNL_CantOpen:
            SCIPerrorMessage("%s\n", writerresult.second.c_str());
            return SCIP_FILECREATEERROR;
         case NLW2_WriteNL_Failed:
            SCIPerrorMessage("%s\n", writerresult.second.c_str());
            return SCIP_WRITEERROR;
         case NLW2_WriteNL_Unset:
         default:
            SCIPerrorMessage("%s\n", writerresult.second.c_str());
            return SCIP_ERROR;
      }

      /* copy temporary file into file */
      strcat(tempname, ".nl");
      tempfile = fopen(tempname, "rb");
      if( tempfile == NULL )
      {
         SCIPerrorMessage("Cannot open temporary file <%s> for reading: error %d\n", tempname, errno);
         return SCIP_FILECREATEERROR;
      }

      char buf[1024];
      int n;
      while( (n=fread(buf, 1, sizeof(buf), tempfile)) != 0 )
         fwrite(buf, 1, n, file);

      fclose(tempfile);
      remove(tempname);
   }
   catch( const mp::UnsupportedError& e )
   {
      SCIPerrorMessage("constraint not writable as AMPL .nl: %s\n", e.what());
      return SCIP_WRITEERROR;
   }
   catch( const mp::Error& e )
   {
      // some other error from ampl/mp
      SCIPerrorMessage("%s\n", e.what());
      return SCIP_WRITEERROR;
   }
   catch( const fmt::SystemError& e )
   {
      // probably a file open error
      SCIPerrorMessage("%s\n", e.what());
      return SCIP_FILECREATEERROR;
   }
   catch( const std::bad_alloc& e )
   {
      SCIPerrorMessage("Out of memory: %s\n", e.what());
      return SCIP_NOMEMORY;
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
   SCIP_CALL( SCIPsetReaderRead(scip, reader, readerReadNl) );
   SCIP_CALL( SCIPsetReaderWrite(scip, reader, readerWriteNl) );

   SCIP_CALL( SCIPincludeExternalCodeInformation(scip, "AMPL/MP 4.0.0", "AMPL .nl file reader library (github.com/ampl/mp)") );

   return SCIP_OKAY;
}

/** writes AMPL solution file
 *
 * problem must have been read with .nl reader
 */
SCIP_RETCODE SCIPwriteSolutionNl(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_PROBDATA* probdata;

   assert(scip != NULL);

   probdata = SCIPgetProbData(scip);
   if( probdata == NULL )
   {
      SCIPerrorMessage("No AMPL nl file read. Cannot write AMPL solution.\n");
      return SCIP_ERROR;
   }

   probdata->filenamestub[probdata->filenamestublen] = '.';
   probdata->filenamestub[probdata->filenamestublen+1] = 's';
   probdata->filenamestub[probdata->filenamestublen+2] = 'o';
   probdata->filenamestub[probdata->filenamestublen+3] = 'l';
   probdata->filenamestub[probdata->filenamestublen+4] = '\0';

   FILE* solfile = fopen(probdata->filenamestub, "w");
   if( solfile == NULL )
   {
      SCIPerrorMessage("could not open file <%s> for writing\n", probdata->filenamestub);
      probdata->filenamestub[probdata->filenamestublen] = '\0';

      return SCIP_WRITEERROR;
   }
   probdata->filenamestub[probdata->filenamestublen] = '\0';

   // see ampl/mp:sol.h:WriteSolFile() (seems buggy, https://github.com/ampl/mp/issues/135) and asl/writesol.c for solution file format
   SCIP_CALL( SCIPprintStatus(scip, solfile) );
   SCIPinfoMessage(scip, solfile, "\n\n");

   SCIPinfoMessage(scip, solfile, "Options\n%d\n", probdata->namplopts);
   for( int i = 0; i < probdata->namplopts; ++i )
      SCIPinfoMessage(scip, solfile, "%d\n", probdata->amplopts[i]);

   bool haveprimal = SCIPgetBestSol(scip) != NULL;
   bool havedual = probdata->islp && SCIPgetStage(scip) == SCIP_STAGE_SOLVED && !SCIPhasPerformedPresolve(scip);

   SCIPinfoMessage(scip, solfile, "%d\n%d\n", probdata->nconss, havedual ? probdata->nconss : 0);
   SCIPinfoMessage(scip, solfile, "%d\n%d\n", probdata->nvars, haveprimal ? probdata->nvars : 0);

   SCIPdebug( SCIPprintSol(scip, SCIPgetBestSol(scip), NULL, TRUE); )

   if( havedual )
      for( int c = 0; c < probdata->nconss; ++c )
      {
         SCIP_CONS* transcons;
         SCIP_Real dualval;

         /* dual solution is created by LP solver and therefore only available for linear constraints */
         SCIP_CALL( SCIPgetTransformedCons(scip, probdata->conss[c], &transcons) );
         assert(transcons == NULL || strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(transcons)), "linear") == 0);

         if( transcons == NULL )
            dualval = 0.0;
         else if( SCIPgetObjsense(scip) == SCIP_OBJSENSE_MINIMIZE )
            dualval = SCIPgetDualsolLinear(scip, transcons);
         else
            dualval = -SCIPgetDualsolLinear(scip, transcons);
         assert(dualval != SCIP_INVALID);

         SCIPinfoMessage(scip, solfile, "%.17g\n", dualval);
      }

   if( haveprimal )
      for( int i = 0; i < probdata->nvars; ++i )
         SCIPinfoMessage(scip, solfile, "%.17g\n", SCIPgetSolVal(scip, SCIPgetBestSol(scip), probdata->vars[i]));

   /* AMPL solve status codes are at https://mp.ampl.com/details.html#_CPPv4N2mp3sol6StatusE
    * (mp::sol::Status enum in amplmp/include/mp/common.h)
    */
   int solve_result_num = mp::sol::FAILURE;
   switch( SCIPgetStatus(scip) )
   {
      case SCIP_STATUS_UNKNOWN:
         break;
      case SCIP_STATUS_USERINTERRUPT:
      case SCIP_STATUS_TERMINATE:
         if( haveprimal )
            solve_result_num = mp::sol::LIMIT_FEAS_INTERRUPT;
         else
            solve_result_num = mp::sol::LIMIT_NO_FEAS_INTERRUPT;
         break;
      case SCIP_STATUS_NODELIMIT:
      case SCIP_STATUS_TOTALNODELIMIT:
      case SCIP_STATUS_STALLNODELIMIT:
         if( haveprimal )
            solve_result_num = mp::sol::LIMIT_FEAS_NODES;
         else
            solve_result_num = mp::sol::LIMIT_NO_FEAS_NODES;
         break;
      case SCIP_STATUS_TIMELIMIT:
         if( haveprimal )
            solve_result_num = mp::sol::LIMIT_FEAS_TIME;
         else
            solve_result_num = mp::sol::LIMIT_NO_FEAS_TIME;
         break;
      case SCIP_STATUS_MEMLIMIT:
         if( haveprimal )
            solve_result_num = mp::sol::LIMIT_FEAS_SOFTMEM;
         else
            solve_result_num = mp::sol::LIMIT_NO_FEAS_SOFTMEM;
         break;
      case SCIP_STATUS_GAPLIMIT:
         /* there is no enum value for gaplimit, so use "work limit" */
         if( haveprimal )
            solve_result_num = mp::sol::LIMIT_FEAS_WORK;
         else
            solve_result_num = mp::sol::LIMIT_NO_FEAS_WORK;
         break;
      case SCIP_STATUS_PRIMALLIMIT:
         solve_result_num = mp::sol::LIMIT_FEAS_BESTOBJ;
         break;
      case SCIP_STATUS_DUALLIMIT:
         if( haveprimal )
            solve_result_num = mp::sol::LIMIT_FEAS_BESTBND;
         else
            solve_result_num = mp::sol::LIMIT_NO_FEAS_BESTBND;
         break;
      case SCIP_STATUS_SOLLIMIT:
         if( haveprimal )
            solve_result_num = mp::sol::LIMIT_FEAS_NUMSOLS;
         else  /* reach solution limit without solution? */
            solve_result_num = mp::sol::LIMIT_NO_FEAS;
         break;
      case SCIP_STATUS_BESTSOLLIMIT:
      case SCIP_STATUS_RESTARTLIMIT:
         /* rare SCIP specific limits that don't map to an AMPL status */
         if( haveprimal )
            solve_result_num = mp::sol::LIMIT_FEAS;
         else
            solve_result_num = mp::sol::LIMIT_NO_FEAS;
         break;
      case SCIP_STATUS_OPTIMAL:
         solve_result_num = mp::sol::SOLVED;
         break;
      case SCIP_STATUS_INFEASIBLE:
         solve_result_num = mp::sol::INFEASIBLE;
         break;
      case SCIP_STATUS_UNBOUNDED:
         if( haveprimal )
            solve_result_num = mp::sol::UNBOUNDED_FEAS;
         else
            solve_result_num = mp::sol::UNBOUNDED_NO_FEAS;
         break;
      case SCIP_STATUS_INFORUNBD:
         solve_result_num = mp::sol::LIMIT_INF_UNB;
         break;
   }
   SCIPinfoMessage(scip, solfile, "objno 0 %d\n", solve_result_num);

   if( fclose(solfile) != 0 )
   {
      SCIPerrorMessage("could not close solution file after writing\n");
      return SCIP_WRITEERROR;
   }

   return SCIP_OKAY;
}
