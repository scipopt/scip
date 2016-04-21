/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2016 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   cons_expr.c
 * @brief  constraint handler for expression constraints (in particular, nonlinear constraints)
 * @author Stefan Vigerske
 * @author Benjamin Mueller
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>
#include <ctype.h>

#include "scip/cons_expr.h"
#include "scip/struct_cons_expr.h"
#include "scip/cons_expr_var.h"
#include "scip/cons_expr_value.h"
#include "scip/cons_expr_sumprod.h"

/* fundamental constraint handler properties */
#define CONSHDLR_NAME          "expr"
#define CONSHDLR_DESC          "constraint handler for expressions"
#define CONSHDLR_ENFOPRIORITY       -60 /**< priority of the constraint handler for constraint enforcing */
#define CONSHDLR_CHECKPRIORITY -4000010 /**< priority of the constraint handler for checking feasibility */
#define CONSHDLR_EAGERFREQ          100 /**< frequency for using all instead of only the useful constraints in separation,
                                         *   propagation and enforcement, -1 for no eager evaluations, 0 for first only */
#define CONSHDLR_NEEDSCONS         TRUE /**< should the constraint handler be skipped, if no constraints are available? */

/* optional constraint handler properties */
#define CONSHDLR_SEPAPRIORITY        10 /**< priority of the constraint handler for separation */
#define CONSHDLR_SEPAFREQ             1 /**< frequency for separating cuts; zero means to separate only in the root node */
#define CONSHDLR_DELAYSEPA        FALSE /**< should separation method be delayed, if other separators found cuts? */

#define CONSHDLR_PROPFREQ             1 /**< frequency for propagating domains; zero means only preprocessing propagation */
#define CONSHDLR_DELAYPROP        FALSE /**< should propagation method be delayed, if other propagators found reductions? */
#define CONSHDLR_PROP_TIMING     SCIP_PROPTIMING_BEFORELP /**< propagation timing mask of the constraint handler*/

#define CONSHDLR_PRESOLTIMING    SCIP_PRESOLTIMING_ALWAYS /**< presolving timing of the constraint handler (fast, medium, or exhaustive) */
#define CONSHDLR_MAXPREROUNDS        -1 /**< maximal number of presolving rounds the constraint handler participates in (-1: no limit) */



/* enable nonlinear constraint upgrading */
#include "scip/cons_nonlinear.h"
#define NONLINCONSUPGD_PRIORITY   100000 /**< priority of the constraint handler for upgrading of nonlinear constraints */



/** ensures that a block memory array has at least a given size
 *
 *  if cursize is 0, then *array1 can be NULL
 */
#define ENSUREBLOCKMEMORYARRAYSIZE(scip, array1, cursize, minsize)      \
   do {                                                                 \
      int __newsize;                                                    \
      assert((scip)  != NULL);                                          \
      if( (cursize) >= (minsize) )                                      \
         break;                                                         \
      __newsize = SCIPcalcMemGrowSize(scip, minsize);                   \
      assert(__newsize >= (minsize));                                   \
      SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &(array1), cursize, __newsize) ); \
      (cursize) = __newsize;                                            \
   } while( FALSE )


/*
 * Data structures
 */

/* TODO: fill in the necessary constraint data */

/** constraint data for expr constraints */
struct SCIP_ConsData
{
   SCIP_CONSEXPR_EXPR*   expr;               /**< expression that represents this constraint (must evaluate to 0 (FALSE) or 1 (TRUE)) */
   SCIP_Real             lhs;                /**< left-hand side */
   SCIP_Real             rhs;                /**< right-hand side */
};

/** constraint handler data */
struct SCIP_ConshdlrData
{
   SCIP_CONSEXPR_EXPRHDLR** exprhdlrs;       /**< expression handlers */
   int                      nexprhdlrs;      /**< number of expression handlers */
   int                      exprhdlrssize;   /**< size of exprhdlrs array */

   SCIP_CONSEXPR_EXPRHDLR*  exprvarhdlr;     /**< variable expression handler */
   SCIP_CONSEXPR_EXPRHDLR*  exprvalhdlr;     /**< value expression handler */
   SCIP_CONSEXPR_EXPRHDLR*  exprsumhdlr;     /**< summation expression handler */
   SCIP_CONSEXPR_EXPRHDLR*  exprprodhdlr;    /**< product expression handler */
};

/* data passed on during expression evaluation in a point */
typedef struct
{
   SCIP_SOL*             sol;                /**< solution that is evaluated */
   unsigned int          soltag;             /**< solution tag */
   SCIP_Bool             aborted;            /**< whether the evaluation has been aborted due to an evaluation error */
} EXPREVAL_DATA;

/* data passed on during expression evaluation over a box */
typedef struct
{
   unsigned int          boxtag;             /**< box tag */
   SCIP_Bool             aborted;            /**< whether the evaluation has been aborted due to an empty interval */
} EXPRINTEVAL_DATA;

/* data passed on during variable locking  */
typedef struct
{
   SCIP_CONSEXPR_EXPRHDLR* exprvarhdlr;      /**< handler for variable expressions (to recognize variable expressions) */
   int                     nlockspos;        /**< number of positive locks */
   int                     nlocksneg;        /**< number of negative locks */
} EXPRLOCK_DATA;

/*
 * Local methods
 */

/*
 * Walking methods: several operations need to traverse the whole expression tree: print, evaluate, free, etc.
 * These operations have a very natural recursive implementation. However, deep recursion can raise stack overflows.
 * To avoid this issue, the method SCIPwalkConsExprExprDF is introduce to traverse the tree and execute callbacks
 * at different places. Here are the callbacks needed for performing the mentioned operations
 */

/** expression walk callback to free an expression including its children (if not used anywhere else) */
static
SCIP_DECL_CONSEXPREXPRWALK_VISIT(freeExpr)
{
   assert(expr != NULL);

   /* the expression should not be used by more than the parent that is also freed */
   assert(expr->nuses == 1);

   switch( stage )
   {
      case SCIP_CONSEXPREXPRWALK_ENTEREXPR :
      {
         /* free expression data, if any, when entering expression */

         if( expr->exprdata != NULL && expr->exprhdlr->freedata != NULL )
         {
            SCIP_CALL( expr->exprhdlr->freedata(scip, expr) );
            assert(expr->exprdata == NULL);
         }

         *result = SCIP_CONSEXPREXPRWALK_CONTINUE;
         return SCIP_OKAY;
      }

      case SCIP_CONSEXPREXPRWALK_VISITINGCHILD :
      {
         /* check whether a child needs to be visited (nuses == 1)
          * if not, then we still have to release it
          */
         SCIP_CONSEXPR_EXPR* child;

         assert(expr->walkcurrentchild < expr->nchildren);

         child = expr->children[expr->walkcurrentchild];
         if( child->nuses > 1 )
         {
            /* child is not going to be freed: just release it */
            SCIP_CALL( SCIPreleaseConsExprExpr(scip, &child) );
            *result = SCIP_CONSEXPREXPRWALK_SKIP;
         }
         else
         {
            assert(child->nuses == 1);
            *result = SCIP_CONSEXPREXPRWALK_CONTINUE;
         }

         return SCIP_OKAY;
      }

      case SCIP_CONSEXPREXPRWALK_LEAVEEXPR :
      {
         /* free expression when leaving it */

         /* free children array, if any */
         SCIPfreeBlockMemoryArrayNull(scip, &expr->children, expr->childrensize);

         SCIPfreeBlockMemory(scip, &expr);
         assert(expr == NULL);

         *result = SCIP_CONSEXPREXPRWALK_CONTINUE;
         return SCIP_OKAY;
      }

      case SCIP_CONSEXPREXPRWALK_VISITEDCHILD :
      default:
      {
         SCIPABORT(); /* we should never be called in this stage */
         *result = SCIP_CONSEXPREXPRWALK_CONTINUE;
         return SCIP_OKAY;
      }
   }
}

/** expression walk callback to print an expression */
static
SCIP_DECL_CONSEXPREXPRWALK_VISIT(printExpr)
{
   FILE* file;

   assert(expr != NULL);
   assert(expr->exprhdlr != NULL);

   file = (FILE*)data;

   if( expr->exprhdlr->print == NULL )
   {
      /* default: <hdlrname>(<child1>, <child2>, ...) */
      switch( stage )
      {
         case SCIP_CONSEXPREXPRWALK_ENTEREXPR :
         {
            SCIPinfoMessage(scip, file, SCIPgetConsExprExprHdlrName(expr->exprhdlr));
            if( expr->nchildren > 0 )
            {
               SCIPinfoMessage(scip, file, "(");
            }
            break;
         }

         case SCIP_CONSEXPREXPRWALK_VISITEDCHILD :
         {
            if( SCIPgetConsExprExprWalkCurrentChild(expr) < expr->nchildren-1 )
            {
               SCIPinfoMessage(scip, file, ", ");
            }
            else
            {
               SCIPinfoMessage(scip, file, ")");
            }

            break;
         }

         default: ;
      }
   }
   else
   {
      /* redirect to expression callback */
      SCIP_CALL( (*expr->exprhdlr->print)(scip, expr, stage, file) );
   }

   *result = SCIP_CONSEXPREXPRWALK_CONTINUE;

   return SCIP_OKAY;
}

/** expression walk callback when evaluating expression, called before child is visited */
static
SCIP_DECL_CONSEXPREXPRWALK_VISIT(evalExprVisitChild)
{
   EXPREVAL_DATA* evaldata;

   assert(expr != NULL);
   assert(data != NULL);

   evaldata = (EXPREVAL_DATA*)data;

   /* skip child if it has been evaluated for that solution already */
   if( evaldata->soltag != 0 && evaldata->soltag == expr->children[expr->walkcurrentchild]->evaltag )
   {
      if( expr->children[expr->walkcurrentchild]->evalvalue == SCIP_INVALID )
      {
         evaldata->aborted = TRUE;
         *result = SCIP_CONSEXPREXPRWALK_ABORT;
      }
      else
      {
         *result = SCIP_CONSEXPREXPRWALK_SKIP;
      }
   }
   else
   {
      *result = SCIP_CONSEXPREXPRWALK_CONTINUE;
   }

   return SCIP_OKAY;
}

/** expression walk callback when evaluating expression, called when expression is left */
static
SCIP_DECL_CONSEXPREXPRWALK_VISIT(evalExprLeaveExpr)
{
   EXPREVAL_DATA* evaldata;

   assert(expr != NULL);
   assert(data != NULL);
   assert(expr->exprhdlr->eval != NULL);

   evaldata = (EXPREVAL_DATA*)data;

   SCIP_CALL( (*expr->exprhdlr->eval)(scip, expr, &expr->evalvalue, evaldata->sol) );
   expr->evaltag = evaldata->soltag;

   if( expr->evalvalue == SCIP_INVALID )
   {
      evaldata->aborted = TRUE;
      *result = SCIP_CONSEXPREXPRWALK_ABORT;
   }
   else
   {
      *result = SCIP_CONSEXPREXPRWALK_CONTINUE;
   }

   return SCIP_OKAY;
}

/** expression walk callback when evaluating expression on intervals, called before child is visited */
static
SCIP_DECL_CONSEXPREXPRWALK_VISIT(intevalExprVisitChild)
{
   EXPRINTEVAL_DATA* propdata;

   assert(expr != NULL);
   assert(data != NULL);

   propdata = (EXPRINTEVAL_DATA*)data;

   /* skip child if it has been evaluated already */
   if( propdata->boxtag != 0 && propdata->boxtag == expr->children[expr->walkcurrentchild]->intevaltag )
   {
      if( SCIPintervalIsEmpty(SCIPinfinity(scip), expr->children[expr->walkcurrentchild]->interval) )
      {
         propdata->aborted = TRUE;
         *result = SCIP_CONSEXPREXPRWALK_ABORT;
      }
      else
      {
         *result = SCIP_CONSEXPREXPRWALK_SKIP;
      }
   }
   else
   {
      *result = SCIP_CONSEXPREXPRWALK_CONTINUE;
   }

   return SCIP_OKAY;
}

/** expression walk callback when evaluating expression on intervals, called when expression is left */
static
SCIP_DECL_CONSEXPREXPRWALK_VISIT(intevalExprLeaveExpr)
{
   EXPRINTEVAL_DATA* propdata;

   assert(expr != NULL);
   assert(data != NULL);

   propdata = (EXPRINTEVAL_DATA*)data;

   /* set tag in any case */
   expr->intevaltag = propdata->boxtag;

   /* set interval to [-inf,+inf] if interval evaluation callback is not implemented */
   if( expr->exprhdlr->inteval == NULL )
   {
      SCIPintervalSetEntire(SCIPinfinity(scip), &expr->interval);
      *result = SCIP_CONSEXPREXPRWALK_CONTINUE;

      return SCIP_OKAY;
   }

   /* evaluate current expression and move on */
   SCIP_CALL( (*expr->exprhdlr->inteval)(scip, expr, &expr->interval) );

   /* stop if callback returned an empty interval */
   if( SCIPintervalIsEmpty(SCIPinfinity(scip), expr->interval) )
   {
      propdata->aborted = TRUE;
      *result = SCIP_CONSEXPREXPRWALK_ABORT;
   }
   else
   {
      *result = SCIP_CONSEXPREXPRWALK_CONTINUE;
   }

   return SCIP_OKAY;
}

static
SCIP_DECL_CONSEXPREXPRWALK_VISIT(lockVar)
{
   EXPRLOCK_DATA* lockdata;

   assert(expr != NULL);
   assert(data != NULL);
   assert(result != NULL);
   assert(stage == SCIP_CONSEXPREXPRWALK_ENTEREXPR);

   lockdata = (EXPRLOCK_DATA*)data;

   if( SCIPgetConsExprExprHdlr(expr) == lockdata->exprvarhdlr )
   {
      /* if a variable, lock in both directions */
      SCIP_CALL( SCIPaddVarLocks(scip, SCIPgetConsExprExprVarVar(expr), lockdata->nlockspos, lockdata->nlocksneg) );
   }

   *result = SCIP_CONSEXPREXPRWALK_CONTINUE;

   return SCIP_OKAY;
}


/** @name Parsing methods
 * @{
 * Here is an attempt at defining the grammar of an expression.
 * We use upper case names for variables (in the grammar sense) and terminals are between "".
 * Loosely speaking, a Base will be any "block", a Factor is a Base to a power, a Term is a product of Factors
 * and an Expression is a sum of terms.
 * The actual definition:
 * <pre>
 * Expression -> ["+" | "-"] Term { ("+" | "-" | "number *") ] Term }
 * Term       -> Factor { ("*" | "/" ) Factor }
 * Factor     -> Base [ "^" "number" | "^(" "number" ")" ]
 * Base       -> "number" | "<varname>" | "(" Expression ")" | Op "(" OpExpression ")
 * </pre>
 * where [a|b] means a or b or none, (a|b) means a or b, {a} means 0 or more a.
 *
 * Note that Op and OpExpression are undefined. Op corresponds to the name of an expression handler and
 * OpExpression to whatever string the expression handler accepts (through its parse method).
 *
 * Each parse(Expr|Factor|Term|Base) returns an SCIP_CONSEXPR_EXPR.
 * For the time being we assume that
 * - Expr is always a exprsum,
 * - Term is always a exprprod,
 * - Factor is always a exprprod,
 * - Base can be anything.
 *
 * @todo We can change the grammar so that Factor becomes base and we allow a Term to be
 *       <pre> Term       -> Factor { ("*" | "/" | "^") Factor } </pre>
 */

#ifdef PARSE_DEBUG
#define debugParse                      printf
#else
#define debugParse                      while( FALSE ) printf
#endif
static
SCIP_RETCODE parseExpr(SCIP*, SCIP_CONSHDLR*, SCIP_HASHMAP*, const char*, const char**, SCIP_CONSEXPR_EXPR**);

/** Parses base to build a value, variable, sum, or function-like ("func(...)") expression.
 * <pre>
 * Base       -> "number" | "<varname>" | "(" Expression ")" | Op "(" OpExpression ")
 * </pre>
 */
static
SCIP_RETCODE parseBase(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< expression constraint handler */
   SCIP_HASHMAP*         vartoexprvarmap,    /**< hashmap to map between SCIP vars and var expressions */
   const char*           expr,               /**< expr that we are parsing */
   const char**          newpos,             /**< buffer to store the position of expr where we finished reading */
   SCIP_CONSEXPR_EXPR**  basetree            /**< buffer to store the expr parsed by Base */
   )
{
   SCIP_VAR* var;

   debugParse("parsing base from %s\n", expr);

   /* ignore whitespace */
   while( isspace((unsigned char)*expr) )
      ++expr;

   if( *expr == '\0' )
   {
      SCIPerrorMessage("Unexpected end of expression string\n");
      return SCIP_READERROR;
   }

   if( *expr == '<' )
   {
      /* parse a variable */
      SCIP_CALL( SCIPparseVarName(scip, expr, &var, (char**)newpos) );

      if( var == NULL )
      {
         SCIPerrorMessage("Could not find variable with name '%s'\n", expr);
         return SCIP_READERROR;
      }
      expr = *newpos;

      /* check if we have already created an expression out of this var */
      if( SCIPhashmapExists(vartoexprvarmap, (void *)var) )
      {
         debugParse("Variable %s has been parsed, capturing its expression\n", SCIPvarGetName(var));
         *basetree = (SCIP_CONSEXPR_EXPR*)SCIPhashmapGetImage(vartoexprvarmap, (void *)var);
         SCIPcaptureConsExprExpr(*basetree);
      }
      else
      {
         debugParse("First time parsing variable %s, creating varexpr and adding it to hashmap\n", SCIPvarGetName(var));
         SCIP_CALL( SCIPcreateConsExprExprVar(scip, conshdlr, basetree, var) );
         SCIP_CALL( SCIPhashmapInsert(vartoexprvarmap, (void*)var, (void*)(*basetree)) );
      }
   }
   else if( *expr == '(' )
   {
      /* parse expression */
      SCIP_CALL( parseExpr(scip, conshdlr, vartoexprvarmap, ++expr, newpos, basetree) );
      expr = *newpos;

      /* expect ')' */
      if( *expr != ')' )
      {
         SCIPerrorMessage("Expected ')', but got <%c> from <%s>\n", *expr, expr);
         SCIP_CALL( SCIPreleaseConsExprExpr(scip, basetree) );
         return SCIP_READERROR;
      }
      ++expr;
      debugParse("Done parsing expression, continue with <%s>\n", expr);
   }
   else if( isdigit(*expr) )
   {
      /* parse number */
      SCIP_Real value;
      if( !SCIPstrToRealValue(expr, &value, (char**)&expr) )
      {
         SCIPerrorMessage("error parsing number from <%s>\n", expr);
         return SCIP_READERROR;
      }
      debugParse("Parsed value %g, creating a value-expression.\n", value);
      SCIP_CALL( SCIPcreateConsExprExprValue(scip, conshdlr, basetree, value) );
   }
   else if( isalpha(*expr) )
   {
      /* a (function) name is coming, should find exprhandler with such name */
      int i;
      int npar;
      const char* init;
      char operatorname[SCIP_MAXSTRLEN];
      SCIP_CONSEXPR_EXPRHDLR* exprhdlr;
      SCIP_Bool insidevarname;
      SCIP_Bool success;

      /* get name */
      i = 0;
      while( *expr != '(' && !isspace((unsigned char)*expr) && *expr != '\0' )
      {
         operatorname[i] = *expr;
         ++expr;
         ++i;
      }
      operatorname[i] = '\0';

      if( *expr != '(' )
      {
         SCIPerrorMessage("Expected '(' after operator name <%s>, but got %s.\n", operatorname, expr);
         return SCIP_READERROR;
      }

      /* search for expression handler */
      exprhdlr = SCIPfindConsExprExprHdlr(conshdlr, operatorname);

      /* check expression handler exists and has a parsing method */
      if( exprhdlr == NULL )
      {
         SCIPerrorMessage("No expression handler with name <%s> found.\n", operatorname);
         return SCIP_READERROR;
      }
      if( exprhdlr->parse == NULL )
      {
         SCIPerrorMessage("Expression handler <%s> has no parsing method.\n", operatorname);
         return SCIP_READERROR;
      }

      /* in case we would need to extract expression between () */
      ++expr;
      init = expr;
      npar = 0;
      insidevarname = FALSE;
      while( *expr != ')' && npar > 0 )
      {
         /* if we achieve the end of string, probably the expression has unbalanced parenthesis */
         if( *expr == '\0' )
         {
            SCIPerrorMessage("Unexpected end of expression string when parsing '%s(%s'\n", operatorname, init);
            return SCIP_READERROR;
         }

         /* count parenthesis if not inside a variable name */
         if( !insidevarname )
         {
            if( *expr == '(' )
               ++npar;
            else if( *expr == ')' )
               --npar;
         }

         /* check if we are entering or leaving a variable name */
         if( insidevarname && *expr == '>' )
            insidevarname = FALSE;
         if( !insidevarname && *expr == '<' )
            insidevarname = TRUE;

         ++expr;
         assert(npar >= 0);
      }
      assert(npar == 0);
      assert(!insidevarname);
      assert(*expr == ')');

      /* temporarily set null-byte and call exprhdlr parser */
      *(char*)expr = '\0';
      SCIP_CALL( exprhdlr->parse(scip, conshdlr, init, basetree, &success) );

      if( !success )
      {
         SCIPerrorMessage("Error while expression handler <%s> was parsing %s\n", operatorname, init);
         assert(*basetree == NULL);
         return SCIP_READERROR;
      }

      *(char*)expr = ')';
   }
   else
   {
      /* Base -> "number" | "<varname>" | "(" Expression ")" | Op "(" OpExpression ") */
      SCIPerrorMessage("Expected a number, (expression), <varname>, Opname(Opexpr), instead got <%c> from %s\n", *expr, expr);
      return SCIP_READERROR;
   }

   *newpos = expr;

   return SCIP_OKAY;
}

/** Parses a factor and builds a product-expression if there is an exponent, otherwise returns the base expression.
 * <pre>
 * Factor -> Base [ "^" "number" | "^(" "number" ")" ]
 * </pre>
 */
static
SCIP_RETCODE parseFactor(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< expression constraint handler */
   SCIP_HASHMAP*         vartoexprvarmap,    /**< hashmap to map between scip vars and var expressions */
   const char*           expr,               /**< expr that we are parsing */
   const char**          newpos,             /**< buffer to store the position of expr where we finished reading */
   SCIP_Real*            exponent,           /**< buffer to store exponent of Factor */
   SCIP_CONSEXPR_EXPR**  factortree          /**< buffer to store the expr parsed by Factor */
   )
{
   SCIP_CONSEXPR_EXPR*  basetree;

   debugParse("parsing factor from %s\n", expr);

   if( *expr == '\0' )
   {
      SCIPerrorMessage("Unexpected end of expression string.\n");
      return SCIP_READERROR;
   }

   /* parse Base */
   /* ignore whitespace */
   while( isspace((unsigned char)*expr) )
      ++expr;

   SCIP_CALL( parseBase(scip, conshdlr, vartoexprvarmap, expr, newpos, &basetree) );
   *factortree = basetree;
   expr = *newpos;

   /* check if there is an exponent */
   /* ignore whitespace */
   while( isspace((unsigned char)*expr) )
      ++expr;
   if( *expr == '^' )
   {

      ++expr;
      while( isspace((unsigned char)*expr) )
         ++expr;

      if( *expr == '\0' )
      {
         SCIPerrorMessage("Unexpected end of expression string after '^'.\n");
         SCIP_CALL( SCIPreleaseConsExprExpr(scip, &basetree) );
         return SCIP_READERROR;
      }

      if( *expr == '(' )
      {
         ++expr;

         /* it is exponent with parenthesis; expect number possibly starting with + or - */
         if( !SCIPstrToRealValue(expr, exponent, (char**)&expr) )
         {
            SCIPerrorMessage("error parsing number from <%s>\n", expr);
            SCIP_CALL( SCIPreleaseConsExprExpr(scip, &basetree) );
            return SCIP_READERROR;
         }

         /* expect the ')' */
         while( isspace((unsigned char)*expr) )
            ++expr;
         if( *expr != ')' )
         {
            SCIPerrorMessage("error in parsing exponent: expected ')', received <%c> from <%s>\n", *expr,  expr);
            SCIP_CALL( SCIPreleaseConsExprExpr(scip, &basetree) );
            return SCIP_READERROR;
         }
         ++expr;
      }
      else
      {
         /* no parenthesis, we should see just a positive number */

         /* expect a digit */
         if( isdigit(*expr) )
         {
            if( !SCIPstrToRealValue(expr, exponent, (char**)&expr) )
            {
               SCIPerrorMessage("error parsing number from <%s>\n", expr);
               SCIP_CALL( SCIPreleaseConsExprExpr(scip, &basetree) );
               return SCIP_READERROR;
            }
         }
         else
         {
            SCIPerrorMessage("error in parsing exponent, expected a digit, received <%c> from <%s>\n", *expr,  expr);
            SCIP_CALL( SCIPreleaseConsExprExpr(scip, &basetree) );
            return SCIP_READERROR;
         }
      }

      debugParse("parsed the exponent %g\n", *exponent);
   }
   else
   {
      /* there is no explicit exponent */
      *exponent = 1.0;
   }

   *newpos = expr;

   return SCIP_OKAY;
}

/** Parses a term and builds a product-expression, where each factor is a child.
 * <pre>
 * Term -> Factor { ("*" | "/" ) Factor }
 * </pre>
 */
static
SCIP_RETCODE parseTerm(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< expression constraint handler */
   SCIP_HASHMAP*         vartoexprvarmap,    /**< hashmap to map between scip vars and var expressions */
   const char*           expr,               /**< expr that we are parsing */
   const char**          newpos,             /**< buffer to store the position of expr where we finished reading */
   SCIP_CONSEXPR_EXPR**  termtree            /**< buffer to store the expr parsed by Term */
   )
{
   SCIP_CONSEXPR_EXPR* factortree;
   SCIP_Real exponent;

   debugParse("parsing term from %s\n", expr);

   /* parse Factor */
   /* ignore whitespace */
   while( isspace((unsigned char)*expr) )
      ++expr;

   SCIP_CALL( parseFactor(scip, conshdlr, vartoexprvarmap, expr, newpos, &exponent, &factortree) );
   expr = *newpos;

   debugParse("back to parsing Term (we have a Factor with exponent %g), continue parsing from %s\n", exponent, expr);

   /* check if Terms has another Factor incoming */
   while( isspace((unsigned char)*expr) )
      ++expr;
   if( *expr == '*' || *expr == '/' )
   {
      /* initialize termtree as a product expression with a single term, so we can append the extra Factors */
      SCIP_CALL( SCIPcreateConsExprExprProduct(scip, conshdlr, termtree, 1, &factortree, &exponent, 1.0) );

      /* loop: parse Factor, find next symbol */
      do
      {
         SCIP_RETCODE retcode;
         SCIP_Bool isdivision;

         isdivision = (*expr == '/') ? TRUE : FALSE;

         debugParse("while parsing term, read char %c\n", *expr);

         ++expr;
         retcode = parseFactor(scip, conshdlr, vartoexprvarmap, expr, newpos, &exponent, &factortree);

         /* release termtree, if parseFactor fails with a read-error */
         if( retcode == SCIP_READERROR )
         {
            SCIP_CALL( SCIPreleaseConsExprExpr(scip, termtree) );
         }
         SCIP_CALL( retcode );

         /* append newly created factor */
         exponent = isdivision ? -exponent : exponent;
         SCIP_CALL( SCIPappendConsExprExprProductExpr(scip, *termtree, factortree, exponent) );

         /* find next symbol */
         expr = *newpos;
         while( isspace((unsigned char)*expr) )
            ++expr;
      } while( *expr == '*' || *expr == '/' );
   }
   else
   {
      /* Term consists of this unique factor^exponent */
      if( exponent != 1.0 )
      {
         SCIP_CALL( SCIPcreateConsExprExprProduct(scip, conshdlr, termtree, 1, &factortree, &exponent, 1.0) );
      }
      else
         *termtree = factortree;
   }

   *newpos = expr;

   return SCIP_OKAY;
}

/** Parses an expression and builds a sum-expression with children.
 * <pre>
 * Expression -> ["+" | "-"] Term { ("+" | "-" | "number *") ] Term }
 * </pre>
 */
static
SCIP_RETCODE parseExpr(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< expression constraint handler */
   SCIP_HASHMAP*         vartoexprvarmap,    /**< hashmap to map between scip vars and var expressions */
   const char*           expr,               /**< expr that we are parsing */
   const char**          newpos,             /**< buffer to store the position of expr where we finished reading */
   SCIP_CONSEXPR_EXPR**  exprtree            /**< buffer to store the expr parsed by Expr */
   )
{
   SCIP_Real sign;
   SCIP_CONSEXPR_EXPR* termtree;

   debugParse("parsing expression %s\n", expr);

   /* ignore whitespace */
   while( isspace((unsigned char)*expr) )
      ++expr;

   /* if '+' or '-', store it */
   sign = 1.0;
   if( *expr == '+' || *expr == '-' )
   {
      debugParse("while parsing expression, read char %c\n", *expr);
      sign = *expr == '+' ? 1.0 : -1.0;
      ++expr;
   }

   SCIP_CALL( parseTerm(scip, conshdlr, vartoexprvarmap, expr, newpos, &termtree) );
   expr = *newpos;

   debugParse("back to parsing expression (we have the following term), continue parsing from %s\n", expr);

   /* check if Expr has another Term incoming */
   while( isspace((unsigned char)*expr) )
      ++expr;
   if( *expr == '+' || *expr == '-' )
   {
      if( SCIPgetConsExprExprHdlr(termtree) == SCIPgetConsExprExprHdlrValue(conshdlr) )
      {
         /* initialize exprtree as a sum expression with a constant only, so we can append the following terms */
         SCIP_CALL( SCIPcreateConsExprExprSum(scip, conshdlr, exprtree, 0, NULL, NULL, sign * SCIPgetConsExprExprValueValue(termtree)) );
         SCIP_CALL( SCIPreleaseConsExprExpr(scip, &termtree) );
      }
      else
      {
         /* initialize exprtree as a sum expression with a single term, so we can append the following terms */
         SCIP_CALL( SCIPcreateConsExprExprSum(scip, conshdlr, exprtree, 1, &termtree, &sign, 0.0) );
      }

      /* loop: parse Term, find next symbol */
      do
      {
         SCIP_RETCODE retcode;
         SCIP_Real coef;

         /* check if we have a "coef * <term>" */
         if( SCIPstrToRealValue(expr, &coef, (char**)newpos) )
         {
            while( isspace((unsigned char)**newpos) )
               ++(*newpos);

            if( **newpos != '*' )
            {
               /* no '*', so fall back to parsing term after sign */
               coef = (*expr == '+') ? 1.0 : -1.0;
               ++expr;
            }
            else
            {
               /* keep coefficient in coef and continue parsing term after coefficient */
               expr = (*newpos)+1;

               while( isspace((unsigned char)*expr) )
                  ++expr;
            }
         }
         else
         {
            coef = (*expr == '+') ? 1.0 : -1.0;
            ++expr;
         }

         debugParse("while parsing expression, read coefficient %g\n", coef);

         retcode = parseTerm(scip, conshdlr, vartoexprvarmap, expr, newpos, &termtree);

         /* release exprtree if parseTerm fails with an read-error */
         if( retcode == SCIP_READERROR )
         {
            SCIP_CALL( SCIPreleaseConsExprExpr(scip, exprtree) );
         }
         SCIP_CALL( retcode );

         /* append newly created term */
         SCIP_CALL( SCIPappendConsExprExprSumExpr(scip, *exprtree, termtree, coef) );

         /* find next symbol */
         expr = *newpos;
         while( isspace((unsigned char)*expr) )
            ++expr;
      } while( *expr == '+' || *expr == '-' );
   }
   else
   {
      /* Expr consists of this unique ['+' | '-'] Term */
      if( sign  < 0.0 )
      {
         assert(sign == -1.0);
         SCIP_CALL( SCIPcreateConsExprExprSum(scip, conshdlr, exprtree, 1, &termtree, &sign, 0.0) );
      }
      else
         *exprtree = termtree;
   }

   *newpos = expr;

   return SCIP_OKAY;
}

/** @} */

/*
 * Callback methods of constraint handler
 */

/* TODO: Implement all necessary constraint handler methods. The methods with #if 0 ... #else #define ... are optional */

static
SCIP_DECL_NONLINCONSUPGD(nonlinconsUpgdExpr)
{
   SCIP_CONSHDLR* consexprhdlr;
   SCIP_EXPRGRAPH* exprgraph;
   SCIP_EXPRGRAPHNODE* node;
   SCIP_CONSEXPR_EXPR* expr;

   assert(nupgdconss != NULL);
   assert(upgdconss != NULL);

   *nupgdconss = 0;

   exprgraph = SCIPgetExprgraphNonlinear(scip, SCIPconsGetHdlr(cons));
   node = SCIPgetExprgraphNodeNonlinear(scip, cons);

   /* no interest in linear constraints */
   if( node == NULL )
      return SCIP_OKAY;

   consexprhdlr = SCIPfindConshdlr(scip, "expr");
   assert(consexprhdlr != NULL);

   /* try to create a cons_expr expression from an expression graph node */
   SCIP_CALL( SCIPcreateConsExprExpr3(scip, consexprhdlr, &expr, exprgraph, node) );

   /* if that didn't work, then because we do not support a certain expression type yet -> no upgrade */
   if( expr == NULL )
      return SCIP_OKAY;

   if( upgdconsssize < 1 )
   {
      /* request larger upgdconss array */
      *nupgdconss = -1;
      return SCIP_OKAY;
   }

   if( SCIPgetNLinearVarsNonlinear(scip, cons) > 0 )
   {
      /* add linear terms */
      SCIP_CONSEXPR_EXPR* varexpr;
      int i;

      /* ensure expr is a sum expression */
      if( SCIPgetConsExprExprHdlr(expr) != SCIPgetConsExprExprHdlrSum(consexprhdlr) )
      {
         SCIP_CONSEXPR_EXPR* sumexpr;

         SCIP_CALL( SCIPcreateConsExprExprSum(scip, consexprhdlr, &sumexpr, 1, &expr, NULL, 0.0) );
         SCIP_CALL( SCIPreleaseConsExprExpr(scip, &expr) );

         expr = sumexpr;
      }

      for( i = 0; i < SCIPgetNLinearVarsNonlinear(scip, cons); ++i )
      {
         SCIP_CALL( SCIPcreateConsExprExprVar(scip, consexprhdlr, &varexpr, SCIPgetLinearVarsNonlinear(scip, cons)[i]) );
         SCIP_CALL( SCIPappendConsExprExprSumExpr(scip, expr, varexpr, SCIPgetLinearCoefsNonlinear(scip, cons)[i]) );
         SCIP_CALL( SCIPreleaseConsExprExpr(scip, &varexpr) );
      }
   }

   *nupgdconss = 1;
   SCIP_CALL( SCIPcreateConsExpr(scip, upgdconss, SCIPconsGetName(cons), expr, SCIPgetLhsNonlinear(scip, cons), SCIPgetRhsNonlinear(scip, cons),
      SCIPconsIsInitial(cons), SCIPconsIsSeparated(cons), SCIPconsIsEnforced(cons),
      SCIPconsIsChecked(cons), SCIPconsIsPropagated(cons), SCIPconsIsLocal(cons),
      SCIPconsIsModifiable(cons), SCIPconsIsDynamic(cons), SCIPconsIsRemovable(cons), FALSE) );

   SCIP_CALL( SCIPreleaseConsExprExpr(scip, &expr) );

   return SCIP_OKAY;
}

/** copy method for constraint handler plugins (called when SCIP copies plugins) */
#if 0
static
SCIP_DECL_CONSHDLRCOPY(conshdlrCopyExpr)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of expr constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define conshdlrCopyExpr NULL
#endif

/** destructor of constraint handler to free constraint handler data (called when SCIP is exiting) */
static
SCIP_DECL_CONSFREE(consFreeExpr)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONSEXPR_EXPRHDLR* exprhdlr;
   int i;

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   for( i = 0; i < conshdlrdata->nexprhdlrs; ++i )
   {
      exprhdlr = conshdlrdata->exprhdlrs[i];
      assert(exprhdlr != NULL);

      if( exprhdlr->freehdlr != NULL )
      {
         SCIP_CALL( (*exprhdlr->freehdlr)(scip, conshdlr, exprhdlr, &exprhdlr->data) );
      }

      SCIPfreeMemory(scip, &exprhdlr->name);
      SCIPfreeMemoryNull(scip, &exprhdlr->desc);

      SCIPfreeMemory(scip, &exprhdlr);
   }

   SCIPfreeBlockMemoryArray(scip, &conshdlrdata->exprhdlrs, conshdlrdata->exprhdlrssize);

   SCIPfreeMemory(scip, &conshdlrdata);
   SCIPconshdlrSetData(conshdlr, NULL);

   return SCIP_OKAY;
}


/** initialization method of constraint handler (called after problem was transformed) */
#if 0
static
SCIP_DECL_CONSINIT(consInitExpr)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of expr constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consInitExpr NULL
#endif


/** deinitialization method of constraint handler (called before transformed problem is freed) */
#if 0
static
SCIP_DECL_CONSEXIT(consExitExpr)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of expr constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consExitExpr NULL
#endif


/** presolving initialization method of constraint handler (called when presolving is about to begin) */
#if 0
static
SCIP_DECL_CONSINITPRE(consInitpreExpr)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of expr constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consInitpreExpr NULL
#endif


/** presolving deinitialization method of constraint handler (called after presolving has been finished) */
#if 0
static
SCIP_DECL_CONSEXITPRE(consExitpreExpr)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of expr constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consExitpreExpr NULL
#endif


/** solving process initialization method of constraint handler (called when branch and bound process is about to begin) */
#if 0
static
SCIP_DECL_CONSINITSOL(consInitsolExpr)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of expr constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consInitsolExpr NULL
#endif


/** solving process deinitialization method of constraint handler (called before branch and bound process data is freed) */
#if 0
static
SCIP_DECL_CONSEXITSOL(consExitsolExpr)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of expr constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consExitsolExpr NULL
#endif


/** frees specific constraint data */
static
SCIP_DECL_CONSDELETE(consDeleteExpr)
{  /*lint --e{715}*/
   assert(consdata != NULL);
   assert(*consdata != NULL);
   assert((*consdata)->expr != NULL);

   SCIP_CALL( SCIPreleaseConsExprExpr(scip, &(*consdata)->expr) );

   SCIPfreeBlockMemory(scip, consdata);

   return SCIP_OKAY;
}


/** transforms constraint data into data belonging to the transformed problem */
#if 0
static
SCIP_DECL_CONSTRANS(consTransExpr)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of expr constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consTransExpr NULL
#endif


/** LP initialization method of constraint handler (called before the initial LP relaxation at a node is solved) */
#if 0
static
SCIP_DECL_CONSINITLP(consInitlpExpr)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of expr constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consInitlpExpr NULL
#endif


/** separation method of constraint handler for LP solutions */
#if 1
static
SCIP_DECL_CONSSEPALP(consSepalpExpr)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of expr constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consSepalpExpr NULL
#endif


/** separation method of constraint handler for arbitrary primal solutions */
#if 0
static
SCIP_DECL_CONSSEPASOL(consSepasolExpr)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of expr constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consSepasolExpr NULL
#endif


/** constraint enforcing method of constraint handler for LP solutions */
static
SCIP_DECL_CONSENFOLP(consEnfolpExpr)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of expr constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}


/** constraint enforcing method of constraint handler for pseudo solutions */
static
SCIP_DECL_CONSENFOPS(consEnfopsExpr)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of expr constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}


/** feasibility check method of constraint handler for integral solutions */
static
SCIP_DECL_CONSCHECK(consCheckExpr)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of expr constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}


/** domain propagation method of constraint handler */
#if 1
static
SCIP_DECL_CONSPROP(consPropExpr)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of expr constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consPropExpr NULL
#endif


/** presolving method of constraint handler */
#if 0
static
SCIP_DECL_CONSPRESOL(consPresolExpr)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of expr constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consPresolExpr NULL
#endif


/** propagation conflict resolving method of constraint handler */
#if 0
static
SCIP_DECL_CONSRESPROP(consRespropExpr)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of expr constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consRespropExpr NULL
#endif


/** variable rounding lock method of constraint handler */
static
SCIP_DECL_CONSLOCK(consLockExpr)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONSDATA* consdata;
   EXPRLOCK_DATA lockdata;

   assert(conshdlr != NULL);
   assert(cons != NULL);

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   if( consdata->expr == NULL )
      return SCIP_OKAY;

   lockdata.exprvarhdlr = conshdlrdata->exprvarhdlr;
   lockdata.nlockspos = nlockspos;
   lockdata.nlocksneg = nlocksneg;

   SCIP_CALL( SCIPwalkConsExprExprDF(scip, consdata->expr, lockVar, NULL, NULL, NULL, &lockdata) );

   return SCIP_OKAY;
}


/** constraint activation notification method of constraint handler */
#if 0
static
SCIP_DECL_CONSACTIVE(consActiveExpr)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of expr constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consActiveExpr NULL
#endif


/** constraint deactivation notification method of constraint handler */
#if 0
static
SCIP_DECL_CONSDEACTIVE(consDeactiveExpr)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of expr constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consDeactiveExpr NULL
#endif


/** constraint enabling notification method of constraint handler */
#if 0
static
SCIP_DECL_CONSENABLE(consEnableExpr)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of expr constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consEnableExpr NULL
#endif


/** constraint disabling notification method of constraint handler */
#if 0
static
SCIP_DECL_CONSDISABLE(consDisableExpr)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of expr constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consDisableExpr NULL
#endif

/** variable deletion of constraint handler */
#if 0
static
SCIP_DECL_CONSDELVARS(consDelvarsExpr)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of expr constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consDelvarsExpr NULL
#endif


/** constraint display method of constraint handler */
static
SCIP_DECL_CONSPRINT(consPrintExpr)
{  /*lint --e{715}*/

   SCIP_CONSDATA* consdata;

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   /* print left hand side for ranged constraints */
   if( !SCIPisInfinity(scip, -consdata->lhs)
      && !SCIPisInfinity(scip, consdata->rhs)
      && !SCIPisEQ(scip, consdata->lhs, consdata->rhs) )
      SCIPinfoMessage(scip, file, "%.15g <= ", consdata->lhs);

   /* print expression */
   if( consdata->expr != NULL )
   {
      SCIP_CALL( SCIPprintConsExprExpr(scip, consdata->expr, file) );
   }
   else
   {
      SCIPinfoMessage(scip, file, "0");
   }

   /* print right hand side */
   if( SCIPisEQ(scip, consdata->lhs, consdata->rhs) )
      SCIPinfoMessage(scip, file, " == %.15g", consdata->rhs);
   else if( !SCIPisInfinity(scip, consdata->rhs) )
      SCIPinfoMessage(scip, file, " <= %.15g", consdata->rhs);
   else if( !SCIPisInfinity(scip, -consdata->lhs) )
      SCIPinfoMessage(scip, file, " >= %.15g", consdata->lhs);
   else
      SCIPinfoMessage(scip, file, " [free]");

   return SCIP_OKAY;
}


/** constraint copying method of constraint handler */
#if 0
static
SCIP_DECL_CONSCOPY(consCopyExpr)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of expr constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consCopyExpr NULL
#endif


/** constraint parsing method of constraint handler */
static
SCIP_DECL_CONSPARSE(consParseExpr)
{  /*lint --e{715}*/
   SCIP_Real  lhs;
   SCIP_Real  rhs;
   const char* endptr;
   SCIP_CONSEXPR_EXPR* consexprtree;

   SCIPdebugMessage("cons_nonlinear::consparse parsing %s\n",str);

   assert(scip != NULL);
   assert(success != NULL);
   assert(str != NULL);
   assert(name != NULL);
   assert(cons != NULL);

   *success = FALSE;

   /* return if string empty */
   if( !*str )
      return SCIP_OKAY;

   endptr = str;

   /* set left and right hand side to their default values */
   lhs = -SCIPinfinity(scip);
   rhs =  SCIPinfinity(scip);

   /* parse constraint to get lhs, rhs, and expression in between (from cons_linear.c::consparse, but parsing whole string first, then getting expression) */

   /* check for left hand side */
   if( isdigit((unsigned char)str[0]) || ((str[0] == '-' || str[0] == '+') && isdigit((unsigned char)str[1])) )
   {
      /* there is a number coming, maybe it is a left-hand-side */
      if( !SCIPstrToRealValue(str, &lhs, (char**)&endptr) )
      {
         SCIPerrorMessage("error parsing number from <%s>\n", str);
         return SCIP_READERROR;
      }

      /* ignore whitespace */
      while( isspace((unsigned char)*endptr) )
         ++endptr;

      if( endptr[0] != '<' || endptr[1] != '=' )
      {
         /* no '<=' coming, so it was the beginning of the expression and not a left-hand-side */
         lhs = -SCIPinfinity(scip);
      }
      else
      {
         /* it was indeed a left-hand-side, so continue parsing after it */
         str = endptr + 2;

         /* ignore whitespace */
         while( isspace((unsigned char)*str) )
            ++str;
      }
   }

   debugParse("str should start at beginning of expr: %s\n", str);

   /* parse expression: so far we did not allocate memory, so can just return in case of readerror */
   SCIP_CALL( SCIPparseConsExprExpr(scip, conshdlr, str, &str, &consexprtree) );

   /* check for left or right hand side */
   while( isspace((unsigned char)*str) )
      ++str;

   /* check for free constraint */
   if( strncmp(str, "[free]", 6) == 0 )
   {
      if( !SCIPisInfinity(scip, -lhs) )
      {
         SCIPerrorMessage("cannot have left hand side and [free] status \n");
         SCIP_CALL( SCIPreleaseConsExprExpr(scip, &consexprtree) );
         return SCIP_OKAY;
      }
      *success = TRUE;
   }
   else
   {
      switch( *str )
      {
         case '<':
            *success = SCIPstrToRealValue(str+2, &rhs, (char**)&endptr);
            break;
         case '=':
            if( !SCIPisInfinity(scip, -lhs) )
            {
               SCIPerrorMessage("cannot have == on rhs if there was a <= on lhs\n");
               SCIP_CALL( SCIPreleaseConsExprExpr(scip, &consexprtree) );
               return SCIP_OKAY;
            }
            else
            {
               *success = SCIPstrToRealValue(str+2, &rhs, (char**)&endptr);
               lhs = rhs;
            }
            break;
         case '>':
            if( !SCIPisInfinity(scip, -lhs) )
            {
               SCIPerrorMessage("cannot have => on rhs if there was a <= on lhs\n");
               SCIP_CALL( SCIPreleaseConsExprExpr(scip, &consexprtree) );
               return SCIP_OKAY;
            }
            else
            {
               *success = SCIPstrToRealValue(str+2, &lhs, (char**)&endptr);
               break;
            }
         case '\0':
            *success = TRUE;
            break;
         default:
            SCIPerrorMessage("unexpected character %c\n", *str);
            SCIP_CALL( SCIPreleaseConsExprExpr(scip, &consexprtree) );
            return SCIP_OKAY;
      }
   }

   /* create constraint */
   SCIP_CALL( SCIPcreateConsExpr(scip, cons, name,
      consexprtree, lhs, rhs,
      initial, separate, enforce, check, propagate, local, modifiable, dynamic, removable, stickingatnode) );
   assert(*cons != NULL);

   SCIP_CALL( SCIPreleaseConsExprExpr(scip, &consexprtree) );

   debugParse("created expression constraint: <%s>\n", SCIPconsGetName(*cons));

   return SCIP_OKAY;
}


/** constraint method of constraint handler which returns the variables (if possible) */
#if 0
static
SCIP_DECL_CONSGETVARS(consGetVarsExpr)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of expr constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consGetVarsExpr NULL
#endif

/** constraint method of constraint handler which returns the number of variables (if possible) */
#if 0
static
SCIP_DECL_CONSGETNVARS(consGetNVarsExpr)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of expr constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consGetNVarsExpr NULL
#endif

/** constraint handler method to suggest dive bound changes during the generic diving algorithm */
#if 0
static
SCIP_DECL_CONSGETDIVEBDCHGS(consGetDiveBdChgsExpr)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of expr constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consGetDiveBdChgsExpr NULL
#endif



/** creates the handler for an expression handler and includes it into the expression constraint handler */
SCIP_RETCODE SCIPincludeConsExprExprHdlrBasic(
   SCIP*                       scip,         /**< SCIP data structure */
   SCIP_CONSHDLR*              conshdlr,     /**< expression constraint handler */
   SCIP_CONSEXPR_EXPRHDLR**    exprhdlr,     /**< buffer where to store expression handler */
   const char*                 name,         /**< name of expression handler (must not be NULL) */
   const char*                 desc,         /**< description of expression handler (can be NULL) */
   unsigned int                precedence,   /**< precedence of expression operation (used for printing) */
   SCIP_DECL_CONSEXPR_EXPREVAL((*eval)),     /**< point evaluation callback (can not be NULL) */
   SCIP_CONSEXPR_EXPRHDLRDATA* data          /**< data of expression handler (can be NULL) */
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata;

   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);
   assert(name != NULL);
   assert(exprhdlr != NULL);

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   SCIP_CALL( SCIPallocClearMemory(scip, exprhdlr) );

   SCIP_CALL( SCIPduplicateMemoryArray(scip, &(*exprhdlr)->name, name, strlen(name)+1) );
   if( desc != NULL )
   {
      SCIP_CALL( SCIPduplicateMemoryArray(scip, &(*exprhdlr)->desc, desc, strlen(desc)+1) );
   }

   (*exprhdlr)->precedence = precedence;
   (*exprhdlr)->eval = eval;
   (*exprhdlr)->data = data;

   ENSUREBLOCKMEMORYARRAYSIZE(scip, conshdlrdata->exprhdlrs, conshdlrdata->exprhdlrssize, conshdlrdata->nexprhdlrs+1);

   conshdlrdata->exprhdlrs[conshdlrdata->nexprhdlrs] = *exprhdlr;
   ++conshdlrdata->nexprhdlrs;

   return SCIP_OKAY;
}

/** set the expression handler callbacks to copy and free an expression handler */
SCIP_RETCODE SCIPsetConsExprExprHdlrCopyFreeHdlr(
   SCIP*                      scip,              /**< SCIP data structure */
   SCIP_CONSHDLR*             conshdlr,          /**< expression constraint handler */
   SCIP_CONSEXPR_EXPRHDLR*    exprhdlr,          /**< expression handler */
   SCIP_DECL_CONSEXPR_EXPRCOPYHDLR((*copyhdlr)), /**< handler copy callback (can be NULL) */
   SCIP_DECL_CONSEXPR_EXPRFREEHDLR((*freehdlr))  /**< handler free callback (can be NULL) */
)
{
   assert(exprhdlr != NULL);

   exprhdlr->copyhdlr = copyhdlr;
   exprhdlr->freehdlr = freehdlr;

   return SCIP_OKAY;
}

/** set the expression handler callbacks to copy and free expression data */
SCIP_RETCODE SCIPsetConsExprExprHdlrCopyFreeData(
   SCIP*                      scip,              /**< SCIP data structure */
   SCIP_CONSHDLR*             conshdlr,          /**< expression constraint handler */
   SCIP_CONSEXPR_EXPRHDLR*    exprhdlr,          /**< expression handler */
   SCIP_DECL_CONSEXPR_EXPRCOPYDATA((*copydata)), /**< expression data copy callback (can be NULL for expressions without data) */
   SCIP_DECL_CONSEXPR_EXPRFREEDATA((*freedata))  /**< expression data free callback (can be NULL if data does not need to be freed) */
)
{
   assert(exprhdlr != NULL);

   exprhdlr->copydata = copydata;
   exprhdlr->freedata = freedata;

   return SCIP_OKAY;
}

/** set the print callback of an expression handler */
SCIP_RETCODE SCIPsetConsExprExprHdlrPrint(
   SCIP*                      scip,          /**< SCIP data structure */
   SCIP_CONSHDLR*             conshdlr,      /**< expression constraint handler */
   SCIP_CONSEXPR_EXPRHDLR*    exprhdlr,      /**< expression handler */
   SCIP_DECL_CONSEXPR_EXPRPRINT((*print))    /**< print callback (can be NULL) */
)
{
   assert(exprhdlr != NULL);

   exprhdlr->print = print;

   return SCIP_OKAY;
}

/** set the parse callback of an expression handler */
SCIP_RETCODE SCIPsetConsExprExprHdlrParse(
   SCIP*                      scip,          /**< SCIP data structure */
   SCIP_CONSHDLR*             conshdlr,      /**< expression constraint handler */
   SCIP_CONSEXPR_EXPRHDLR*    exprhdlr,      /**< expression handler */
   SCIP_DECL_CONSEXPR_EXPRPARSE((*parse))    /**< parse callback (can be NULL) */
)
{
   assert(exprhdlr != NULL);

   exprhdlr->parse = parse;

   return SCIP_OKAY;
}

/** set the interval evaluation callback of an expression handler */
SCIP_RETCODE SCIPsetConsExprExprHdlrIntEval(
   SCIP*                      scip,          /**< SCIP data structure */
   SCIP_CONSHDLR*             conshdlr,      /**< expression constraint handler */
   SCIP_CONSEXPR_EXPRHDLR*    exprhdlr,      /**< expression handler */
   SCIP_DECL_CONSEXPR_EXPRINTEVAL((*inteval))/**< interval evaluation callback (can be NULL) */
)
{
   assert(exprhdlr != NULL);

   exprhdlr->inteval = inteval;

   return SCIP_OKAY;
}

/** gives expression handlers */
SCIP_CONSEXPR_EXPRHDLR** SCIPgetConsExprExprHdlrs(
   SCIP_CONSHDLR*             conshdlr       /**< expression constraint handler */
)
{
   assert(conshdlr != NULL);

   return SCIPconshdlrGetData(conshdlr)->exprhdlrs;
}

/** gives number of expression handlers */
int SCIPgetConsExprExprNHdlrs(
   SCIP_CONSHDLR*             conshdlr       /**< expression constraint handler */
)
{
   assert(conshdlr != NULL);

   return SCIPconshdlrGetData(conshdlr)->nexprhdlrs;
}

/** returns an expression handler of a given name (or NULL if not found) */
SCIP_CONSEXPR_EXPRHDLR* SCIPfindConsExprExprHdlr(
   SCIP_CONSHDLR*             conshdlr,      /**< expression constraint handler */
   const char*                name           /**< name of expression handler */
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata;
   int h;

   assert(conshdlr != NULL);
   assert(name != NULL);

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   for( h = 0; h < conshdlrdata->nexprhdlrs; ++h )
      if( strcmp(SCIPgetConsExprExprHdlrName(conshdlrdata->exprhdlrs[h]), name) == 0 )
         return conshdlrdata->exprhdlrs[h];

   return NULL;
}

/** returns expression handler for variable expressions */
SCIP_CONSEXPR_EXPRHDLR* SCIPgetConsExprExprHdlrVar(
   SCIP_CONSHDLR*             conshdlr       /**< expression constraint handler */
   )
{
   assert(conshdlr != NULL);

   return SCIPconshdlrGetData(conshdlr)->exprvarhdlr;
}

/** returns expression handler for constant value expressions */
SCIP_CONSEXPR_EXPRHDLR* SCIPgetConsExprExprHdlrValue(
   SCIP_CONSHDLR*             conshdlr       /**< expression constraint handler */
   )
{
   assert(conshdlr != NULL);

   return SCIPconshdlrGetData(conshdlr)->exprvalhdlr;
}

/** returns expression handler for sum expressions */
SCIP_CONSEXPR_EXPRHDLR* SCIPgetConsExprExprHdlrSum(
   SCIP_CONSHDLR*             conshdlr       /**< expression constraint handler */
   )
{
   assert(conshdlr != NULL);

   return SCIPconshdlrGetData(conshdlr)->exprsumhdlr;
}

/** returns expression handler for product expressions */
SCIP_CONSEXPR_EXPRHDLR* SCIPgetConsExprExprHdlrProduct(
   SCIP_CONSHDLR*             conshdlr       /**< expression constraint handler */
   )
{
   assert(conshdlr != NULL);

   return SCIPconshdlrGetData(conshdlr)->exprprodhdlr;
}


/** gives the name of an expression handler */
const char* SCIPgetConsExprExprHdlrName(
   SCIP_CONSEXPR_EXPRHDLR*    exprhdlr       /**< expression handler */
)
{
   assert(exprhdlr != NULL);

   return exprhdlr->name;
}

/** gives the description of an expression handler (can be NULL) */
const char* SCIPgetConsExprExprHdlrDescription(
   SCIP_CONSEXPR_EXPRHDLR*    exprhdlr       /**< expression handler */
)
{
   assert(exprhdlr != NULL);

   return exprhdlr->desc;
}

/** gives the precedence of an expression handler */
unsigned int SCIPgetConsExprExprHdlrPrecedence(
   SCIP_CONSEXPR_EXPRHDLR*    exprhdlr       /**< expression handler */
)
{
   assert(exprhdlr != NULL);

   return exprhdlr->precedence;
}

/** gives the data of an expression handler */
SCIP_CONSEXPR_EXPRHDLRDATA* SCIPgetConsExprExprHdlrData(
   SCIP_CONSEXPR_EXPRHDLR*    exprhdlr      /**< expression handler */
)
{
   assert(exprhdlr != NULL);

   return exprhdlr->data;
}

/** creates and captures an expression with given expression data and children */
SCIP_RETCODE SCIPcreateConsExprExpr(
   SCIP*                   scip,             /**< SCIP data structure */
   SCIP_CONSHDLR*          consexprhdlr,     /**< expression constraint handler */
   SCIP_CONSEXPR_EXPR**    expr,             /**< pointer where to store expression */
   SCIP_CONSEXPR_EXPRHDLR* exprhdlr,         /**< expression handler */
   SCIP_CONSEXPR_EXPRDATA* exprdata,         /**< expression data (expression assumes ownership) */
   int                     nchildren,        /**< number of children */
   SCIP_CONSEXPR_EXPR**    children          /**< children (can be NULL if nchildren is 0) */
   )
{
   int c;

   assert(expr != NULL);
   assert(exprhdlr != NULL);
   assert(children != NULL || nchildren == 0);

   SCIP_CALL( SCIPallocClearBlockMemory(scip, expr) );

   (*expr)->exprhdlr = exprhdlr;
   (*expr)->exprdata = exprdata;

   /* initialize an empty interval for interval evaluation */
   SCIPintervalSetEntire(SCIPinfinity(scip), &(*expr)->interval);

   if( nchildren > 0 )
   {
      SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &(*expr)->children, children, nchildren) );
      (*expr)->nchildren = nchildren;
      (*expr)->childrensize = nchildren;

      for( c = 0; c < nchildren; ++c )
         SCIPcaptureConsExprExpr((*expr)->children[c]);
   }

   SCIPcaptureConsExprExpr(*expr);

   return SCIP_OKAY;
}

/** creates and captures an expression with up to two children */
SCIP_RETCODE SCIPcreateConsExprExpr2(
   SCIP*                   scip,             /**< SCIP data structure */
   SCIP_CONSHDLR*          consexprhdlr,     /**< expression constraint handler */
   SCIP_CONSEXPR_EXPR**    expr,             /**< pointer where to store expression */
   SCIP_CONSEXPR_EXPRHDLR* exprhdlr,         /**< expression handler */
   SCIP_CONSEXPR_EXPRDATA* exprdata,         /**< expression data */
   SCIP_CONSEXPR_EXPR*     child1,           /**< first child (can be NULL) */
   SCIP_CONSEXPR_EXPR*     child2            /**< second child (can be NULL) */
   )
{
   assert(expr != NULL);
   assert(exprhdlr != NULL);
   assert(child1 != NULL);
   assert(child2 != NULL);

   if( child1 != NULL && child2 != NULL )
   {
      SCIP_CONSEXPR_EXPR* pair[2] = {child1, child2};

      SCIP_CALL( SCIPcreateConsExprExpr(scip, consexprhdlr, expr, exprhdlr, exprdata, 2, pair) );
   }
   else if( child2 == NULL )
   {
      SCIP_CALL( SCIPcreateConsExprExpr(scip, consexprhdlr, expr, exprhdlr, exprdata, child1 == NULL ? 0 : 1, &child1) );
   }
   else
   {
      /* child2 != NULL, child1 == NULL */
      SCIP_CALL( SCIPcreateConsExprExpr(scip, consexprhdlr, expr, exprhdlr, exprdata, 1, &child2) );
   }

   return SCIP_OKAY;
}

/** creates and captures an expression from a node in an (old-style) expression graph */
SCIP_RETCODE SCIPcreateConsExprExpr3(
   SCIP*                   scip,             /**< SCIP data structure */
   SCIP_CONSHDLR*          consexprhdlr,     /**< expression constraint handler */
   SCIP_CONSEXPR_EXPR**    expr,             /**< pointer where to store expression */
   SCIP_EXPRGRAPH*         exprgraph,        /**< expression graph */
   SCIP_EXPRGRAPHNODE*     node              /**< expression graph node */
   )
{
   SCIP_CONSEXPR_EXPR** children = NULL;
   int nchildren;
   int c = 0;

   assert(expr != NULL);
   assert(node != NULL);

   *expr = NULL;
   nchildren = SCIPexprgraphGetNodeNChildren(node);

   if( nchildren > 0 )
   {
      SCIP_CALL( SCIPallocBufferArray(scip, &children, nchildren) );

      for( c = 0; c < nchildren; ++c )
      {
         SCIP_CALL( SCIPcreateConsExprExpr3(scip, consexprhdlr, &children[c], exprgraph, SCIPexprgraphGetNodeChildren(node)[c]) );
         if( children[c] == NULL )
            goto TERMINATE;
      }

   }

   switch( SCIPexprgraphGetNodeOperator(node) )
   {
      case SCIP_EXPR_CONST :
         SCIP_CALL( SCIPcreateConsExprExprValue(scip, consexprhdlr, expr, SCIPexprgraphGetNodeOperatorReal(node)) );
         break;

      case SCIP_EXPR_VARIDX :
      {
         int varidx;

         varidx = SCIPexprgraphGetNodeOperatorIndex(node);
         assert(varidx >= 0);
         assert(varidx < SCIPexprgraphGetNVars(exprgraph));

         SCIP_CALL( SCIPcreateConsExprExprVar(scip, consexprhdlr, expr, SCIPexprgraphGetVars(exprgraph)[varidx]) );

         break;
      }

      case SCIP_EXPR_PLUS:
      {
         assert(nchildren == 2);
         SCIP_CALL( SCIPcreateConsExprExprSum(scip, consexprhdlr, expr, 2, children, NULL, 0.0) );

         break;
      }

      case SCIP_EXPR_MINUS:
      {
         SCIP_Real coefs[2] = {1.0, -1.0};

         assert(nchildren == 2);
         SCIP_CALL( SCIPcreateConsExprExprSum(scip, consexprhdlr, expr, 2, children, coefs, 0.0) );

         break;
      }

      case SCIP_EXPR_MUL:
      {
         assert(nchildren == 2);

         SCIP_CALL( SCIPcreateConsExprExprProduct(scip, consexprhdlr, expr, 2, children, NULL, 1.0) );

         break;
      }

      case SCIP_EXPR_DIV:
      {
         SCIP_Real coefs[2] = {1.0, -1.0};

         assert(nchildren == 2);
         SCIP_CALL( SCIPcreateConsExprExprProduct(scip, consexprhdlr, expr, 2, children, coefs, 1.0) );

         break;
      }

      case SCIP_EXPR_SQUARE:
      {
         SCIP_Real two = 2.0;

         assert(nchildren == 1);
         SCIP_CALL( SCIPcreateConsExprExprProduct(scip, consexprhdlr, expr, 1, children, &two, 1.0) );

         break;
      }

      case SCIP_EXPR_SQRT:
      {
         SCIP_Real half = 0.5;

         assert(nchildren == 1);
         SCIP_CALL( SCIPcreateConsExprExprProduct(scip, consexprhdlr, expr, 1, children, &half, 1.0) );

         break;
      }

      case SCIP_EXPR_REALPOWER:
      {
         SCIP_Real exponent;

         exponent = SCIPexprgraphGetNodeOperatorReal(node);

         assert(nchildren == 1);
         SCIP_CALL( SCIPcreateConsExprExprProduct(scip, consexprhdlr, expr, 1, children, &exponent, 1.0) );

         break;
      }

      case SCIP_EXPR_INTPOWER:
      {
         SCIP_Real exponent;

         exponent = (SCIP_Real)SCIPexprgraphGetNodeOperatorIndex(node);

         assert(nchildren == 1);
         SCIP_CALL( SCIPcreateConsExprExprProduct(scip, consexprhdlr, expr, 1, children, &exponent, 1.0) );

         break;
      }

      case SCIP_EXPR_SUM:
      {
         SCIP_CALL( SCIPcreateConsExprExprSum(scip, consexprhdlr, expr, nchildren, children, NULL, 0.0) );

         break;
      }

      case SCIP_EXPR_PRODUCT:
      {
         SCIP_CALL( SCIPcreateConsExprExprProduct(scip, consexprhdlr, expr, nchildren, children, NULL, 0.0) );

         break;
      }

      case SCIP_EXPR_LINEAR:
      {
         SCIP_CALL( SCIPcreateConsExprExprSum(scip, consexprhdlr, expr, nchildren, children, SCIPexprgraphGetNodeLinearCoefs(node), SCIPexprgraphGetNodeLinearConstant(node)) );

         break;
      }

      case SCIP_EXPR_QUADRATIC:
      {
         SCIP_QUADELEM quadelem;
         SCIP_CONSEXPR_EXPR* prod;
         int i;

         SCIP_CALL( SCIPcreateConsExprExprSum(scip, consexprhdlr, expr, 0, NULL, NULL, SCIPexprgraphGetNodeQuadraticConstant(node)) );

         /* append linear terms */
         if( SCIPexprgraphGetNodeQuadraticLinearCoefs(node) != NULL )
         {
            for( i = 0; i < nchildren; ++i )
            {
               if( SCIPexprgraphGetNodeQuadraticLinearCoefs(node)[i] != 0.0 )
               {
                  SCIP_CALL( SCIPappendConsExprExprSumExpr(scip, *expr, children[i], SCIPexprgraphGetNodeQuadraticLinearCoefs(node)[i]) );
               }
            }
         }

         /* append quadratic terms */
         for( i = 0; i < SCIPexprgraphGetNodeQuadraticNQuadElements(node); ++i )
         {
            quadelem = SCIPexprgraphGetNodeQuadraticQuadElements(node)[i];

            if( quadelem.idx1 == quadelem.idx2 )
            {
               SCIP_Real two = 2.0;
               SCIP_CALL( SCIPcreateConsExprExprProduct(scip, consexprhdlr, &prod, 1, &children[quadelem.idx1], &two, 1.0) );
            }
            else
            {
               SCIP_CONSEXPR_EXPR* prodchildren[2] = {children[quadelem.idx1], children[quadelem.idx2]};
               SCIP_CALL( SCIPcreateConsExprExprProduct(scip, consexprhdlr, &prod, 2, prodchildren, NULL, 1.0) );
            }

            SCIP_CALL( SCIPappendConsExprExprSumExpr(scip, *expr, prod, quadelem.coef) );

            SCIP_CALL( SCIPreleaseConsExprExpr(scip, &prod) );
         }

         break;
      }

      case SCIP_EXPR_POLYNOMIAL :
      {
         SCIP_EXPRDATA_MONOMIAL* monom;
         int m;

         SCIP_CALL( SCIPcreateConsExprExprSum(scip, consexprhdlr, expr, 0, NULL, NULL, SCIPexprgraphGetNodePolynomialConstant(node)) );

         /* append monomials */
         for( m = 0; m < SCIPexprgraphGetNodePolynomialNMonomials(node); ++m )
         {
            SCIP_Real* exponents;

            monom = SCIPexprgraphGetNodePolynomialMonomials(node)[m];
            exponents = SCIPexprGetMonomialExponents(monom);

            if( SCIPexprGetMonomialNFactors(monom) == 1 && (exponents == NULL || exponents[0] == 1.0) )
            {
               /* monom is linear in child -> append child itself */
               SCIP_CALL( SCIPappendConsExprExprSumExpr(scip, *expr, children[SCIPexprGetMonomialChildIndices(monom)[0]], SCIPexprGetMonomialCoef(monom)) );
            }
            else
            {
               /* monom is nonlinear -> translate into a product expression */
               SCIP_CONSEXPR_EXPR* monomial;
               int f;

               SCIP_CALL( SCIPcreateConsExprExprProduct(scip, consexprhdlr, &monomial, 0, NULL, NULL, 1.0) );

               for( f = 0; f < SCIPexprGetMonomialNFactors(monom); ++f )
               {
                  SCIP_CALL( SCIPappendConsExprExprProductExpr(scip, monomial, children[SCIPexprGetMonomialChildIndices(monom)[f]], exponents ? exponents[f] : 1.0) );
               }

               SCIP_CALL( SCIPappendConsExprExprSumExpr(scip, *expr, monomial, SCIPexprGetMonomialCoef(monom)) );
               SCIP_CALL( SCIPreleaseConsExprExpr(scip, &monomial) );
            }
         }

         break;
      }

      case SCIP_EXPR_SIGNPOWER:
      case SCIP_EXPR_EXP:
      case SCIP_EXPR_LOG:
      case SCIP_EXPR_SIN:
      case SCIP_EXPR_COS:
      case SCIP_EXPR_TAN:
      case SCIP_EXPR_MIN:
      case SCIP_EXPR_MAX:
      case SCIP_EXPR_ABS:
      case SCIP_EXPR_SIGN:
      case SCIP_EXPR_USER:
      default:
         goto TERMINATE;
   }


TERMINATE:
   /* release all created children expressions (c-1...0) */
   for( --c; c >= 0; --c )
   {
      SCIP_CALL( SCIPreleaseConsExprExpr(scip, &children[c]) );
   }

   SCIPfreeBufferArrayNull(scip, &children);

   return SCIP_OKAY;
}

/** captures an expression (increments usage count) */
void SCIPcaptureConsExprExpr(
   SCIP_CONSEXPR_EXPR*   expr               /**< expression */
   )
{
   assert(expr != NULL);

   ++expr->nuses;
}

/** releases an expression (decrements usage count and possibly frees expression) */
SCIP_RETCODE SCIPreleaseConsExprExpr(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSEXPR_EXPR**  expr                /**< pointer to expression to be released */
   )
{
   assert(expr != NULL);
   assert(*expr != NULL);

   if( (*expr)->nuses == 1 )
   {
      SCIP_CALL( SCIPwalkConsExprExprDF(scip, *expr, freeExpr, freeExpr, NULL, freeExpr, NULL) );
      *expr = NULL;

      return SCIP_OKAY;
   }

   --(*expr)->nuses;
   assert((*expr)->nuses > 0);

   return SCIP_OKAY;
}

/** gives the number of children of an expression */
int SCIPgetConsExprExprNChildren(
   SCIP_CONSEXPR_EXPR*   expr               /**< expression */
   )
{
   assert(expr != NULL);

   return expr->nchildren;
}

/** gives the children of an expression (can be NULL if no children) */
SCIP_CONSEXPR_EXPR** SCIPgetConsExprExprChildren(
   SCIP_CONSEXPR_EXPR*   expr               /**< expression */
   )
{
   assert(expr != NULL);

   return expr->children;
}

/** gets the handler of an expression
 *
 * This identifies the type of the expression (sum, variable, ...).
 */
SCIP_CONSEXPR_EXPRHDLR* SCIPgetConsExprExprHdlr(
   SCIP_CONSEXPR_EXPR*   expr                /**< expression */
   )
{
   assert(expr != NULL);

   return expr->exprhdlr;
}

/** gets the expression data of an expression */
SCIP_CONSEXPR_EXPRDATA* SCIPgetConsExprExprData(
   SCIP_CONSEXPR_EXPR*   expr                /**< expression */
   )
{
   assert(expr != NULL);

   return expr->exprdata;
}

/** sets the expression data of an expression
 *
 * The pointer to possible old data is overwritten and the
 * freedata-callback is not called before.
 * This function is intended to be used by expression handler.
 */
void SCIPsetConsExprExprData(
   SCIP_CONSEXPR_EXPR*     expr,             /**< expression */
   SCIP_CONSEXPR_EXPRDATA* exprdata          /**< expression data to be set (can be NULL) */
   )
{
   assert(expr != NULL);

   expr->exprdata = exprdata;
}

/** print an expression as info-message */
SCIP_RETCODE SCIPprintConsExprExpr(
   SCIP*                   scip,             /**< SCIP data structure */
   SCIP_CONSEXPR_EXPR*     expr,             /**< expression to be printed */
   FILE*                   file              /**< file to print to, or NULL for stdout */
   )
{
   assert(expr != NULL);

   SCIP_CALL( SCIPwalkConsExprExprDF(scip, expr, printExpr, printExpr, printExpr, printExpr, (void*)file) );

   return SCIP_OKAY;
}

/** evaluate an expression in a point
 *
 * Initiates an expression walk to also evaluate children, if necessary.
 * Value can be received via SCIPgetConsExprExprEvalValue().
 * If an evaluation error (division by zero, ...) occurs, this value will
 * be set to SCIP_INVALID.
 *
 * If a nonzero \p soltag is passed, then only (sub)expressions are
 * reevaluated that have a different solution tag. If a soltag of 0
 * is passed, then subexpressions are always reevaluated.
 * The tag is stored together with the value and can be received via
 * SCIPgetConsExprExprEvalTag().
 */
SCIP_RETCODE SCIPevalConsExprExpr(
   SCIP*                   scip,             /**< SCIP data structure */
   SCIP_CONSEXPR_EXPR*     expr,             /**< expression to be evaluated */
   SCIP_SOL*               sol,              /**< solution to be evaluated */
   unsigned int            soltag            /**< tag that uniquely identifies the solution (with its values), or 0. */
   )
{
   EXPREVAL_DATA evaldata;

   /* if value is up-to-date, then nothing to do */
   if( soltag != 0 && expr->evaltag == soltag )
      return SCIP_OKAY;

   evaldata.sol = sol;
   evaldata.soltag = soltag;
   evaldata.aborted = FALSE;

   SCIP_CALL( SCIPwalkConsExprExprDF(scip, expr, NULL, evalExprVisitChild, NULL, evalExprLeaveExpr, &evaldata) );

   if( evaldata.aborted )
   {
      expr->evalvalue = SCIP_INVALID;
      expr->evaltag = soltag;
   }

   return SCIP_OKAY;
}

/** evaluates an expression over a box
 *
 * Initiates an expression walk to also evaluate children, if necessary.
 * The resulting interval can be received via SCIPgetConsExprExprEvalInterval().
 * If the box does not overlap with the domain of the function behind the expression
 * (e.g., sqrt([-2,-1]) or 1/[0,0]) this interval will be empty.
 *
 * For variables, the local variable bounds are used as interval.
 *
 * If a nonzero \p boxtag is passed, then only (sub)expressions are
 * reevaluated that have a different tag. If a tag of 0 is passed,
 * then subexpressions are always reevaluated.
 * The tag is stored together with the interval and can be received via
 * SCIPgetConsExprExprEvalIntervalTag().
 */
SCIP_RETCODE SCIPevalConsExprExprInterval(
   SCIP*                   scip,             /**< SCIP data structure */
   SCIP_CONSEXPR_EXPR*     expr,             /**< expression to be evaluated */
   unsigned int            boxtag            /**< tag that uniquely identifies the current variable domains (with its values), or 0 */
   )
{
   EXPRINTEVAL_DATA propdata;

   /* if value is up-to-date, then nothing to do */
   if( boxtag != 0 && expr->intevaltag == boxtag )
      return SCIP_OKAY;

   propdata.aborted = FALSE;
   propdata.boxtag = boxtag;

   SCIP_CALL( SCIPwalkConsExprExprDF(scip, expr, NULL, intevalExprVisitChild, NULL, intevalExprLeaveExpr, &propdata) );

   if( propdata.aborted )
   {
      SCIPintervalSetEmpty(&expr->interval);
      expr->intevaltag = boxtag;
   }

   return SCIP_OKAY;
}

/** gives the value from the last evaluation of an expression (or SCIP_INVALID if there was an eval error) */
SCIP_Real SCIPgetConsExprExprValue(
   SCIP_CONSEXPR_EXPR*     expr              /**< expression */
   )
{
   assert(expr != NULL);

   return expr->evalvalue;
}

/** returns the interval from the last interval evaluation of an expression (interval can be empty) */
SCIP_INTERVAL SCIPgetConsExprExprInterval(
   SCIP_CONSEXPR_EXPR*     expr              /**< expression */
   )
{
   assert(expr != NULL);

   return expr->interval;
}

/** gives the evaluation tag from the last evaluation, or 0 */
unsigned int SCIPgetConsExprExprEvalTag(
   SCIP_CONSEXPR_EXPR*     expr              /**< expression */
   )
{
   assert(expr != NULL);

   return expr->evaltag;
}

/** gives the box tag from the last interval evaluation, or 0 */
unsigned int SCIPgetConsExprExprEvalIntervalTag(
   SCIP_CONSEXPR_EXPR*     expr              /**< expression */
   )
{
   assert(expr != NULL);

   return expr->intevaltag;
}

/** sets the evaluation value */
void SCIPsetConsExprExprEvalValue(
   SCIP_CONSEXPR_EXPR*     expr,             /**< expression */
   SCIP_Real               value,            /**< value to set */
   unsigned int            tag               /**< tag of solution that was evaluated, or 0 */
   )
{
   assert(expr != NULL);

   expr->evalvalue = value;
   expr->evaltag = tag;
}

/** sets the evaluation interval */
void SCIPsetConsExprExprEvalInterval(
   SCIP_CONSEXPR_EXPR*     expr,             /**< expression */
   SCIP_INTERVAL*          interval,         /**< interval to set */
   unsigned int            tag               /**< tag of variable domains that were evaluated, or 0. */
   )
{
   assert(expr != NULL);

   SCIPintervalSetBounds(&expr->interval, SCIPintervalGetInf(*interval), SCIPintervalGetSup(*interval));
   expr->intevaltag = tag;
}

/** walks the expression graph in depth-first manner and executes callbacks at certain places
 *
 * Many algorithms over expression trees need to traverse the tree in depth-first manner and a
 * natural way of implementing this algorithms is using recursion.
 * In general, a function which traverses the tree in depth-first looks like
 * <pre>
 * fun( expr )
 *    enterexpr()
 *    continue skip or abort
 *       for( child in expr->children )
 *          visitingchild()
 *          continue skip or abort
 *          fun(child, data, proceed)
 *          visitedchild()
 *          continue skip or abort
 *    leaveexpr()
 * </pre>
 * Given that some expressions might be quite deep we provide this functionality in an iterative fashion.
 *
 * Consider an expression (x*y) + z + log(x-y).
 * The corresponding expression graph is
 * <pre>
 *           [+]
 *       /    |   \
 *    [*]     |    [log]
 *    / \     |      |
 *   /   \    |     [-]
 *   |   |    |     / \
 *  [x] [y]  [z]  [x] [y]
 * </pre>
 * (where [x] and [y] are actually the same expression).
 *
 * If given a pointer to the [+] expression is given as root to this expression, it will walk
 * the graph in a depth-first manner and call the given callback methods at various stages.
 * - When entering an expression, it calls the enterexpr callback.
 *   The SCIPgetConsExprExprWalkParent() function indicates from where the expression has been entered (NULL for the root expression).
 * - Before visiting a child of an expression, it calls the visitingchild callback.
 *   The SCIPgetConsExprExprWalkCurrentChild() function returns which child will be visited (as an index in the current expr's children array).
 * - When returning from visiting a child of an expression, the visitedchild callback is called.
 *   Again the SCIPgetConsExprExprWalkCurrentChild() function returns which child has been visited.
 * - When leaving an expression, it calls the leaveexpr callback.
 *
 * Thus, for the above expression, the callbacks are called in the following order:
 * - enterexpr([+])
 * - visitingchild([+])  currentchild == 0
 * - enterexpr([*])
 * - visitingchild([*])  currentchild == 0
 * - enterexpr([x])
 * - leaveexpr([x])
 * - visitedchild([*])   currentchild == 0
 * - visitingchild([*])  currentchild == 1
 * - enterexpr([y])
 * - leaveexpr([y])
 * - visitedchild([*])   currentchild == 1
 * - leaveexpr([*])
 * - visitedchild([+])   currentchild == 0
 * - visitingchild([+])  currentchild == 1
 * - enterexpr([z])
 * - leaveexpr([z])
 * - visitedchild([+])   currentchild == 1
 * - visitingchild([+])  currentchild == 2
 * - enterexpr([log])
 * - visitingchild([log]) currentchild == 0
 * - enterexpr([-])
 * - visitingchild([-])  currentchild == 0
 * - enterexpr([x])
 * - leaveexpr([x])
 * - visitedchild([-])   currentchild == 0
 * - visitingchild([-])  currentchild == 1
 * - enterexpr([y])
 * - leaveexpr([y])
 * - visitedchild([-])   currentchild == 1
 * - leaveexpr([-])
 * - visitedchild([log]) currentchild == 0
 * - leaveexpr([log])
 * - visitedchild([+])   currentchild == 2
 * - leaveexpr([+])
 *
 * The callbacks can direct the walking method to skip parts of the tree or abort.
 * If returning SCIP_CONSEXPREXPRWALK_SKIP as result of an enterexpr callback, all children of that expression will be skipped. The leaveexpr callback will still be called.
 * If returning SCIP_CONSEXPREXPRWALK_SKIP as result of an visitingchild callback, visiting the current child will be skipped.
 * If returning SCIP_CONSEXPREXPRWALK_ABORT in any of the callbacks, the walk will be aborted immediately.
 */
SCIP_RETCODE SCIPwalkConsExprExprDF(
   SCIP*                 scip,                         /**< SCIP data structure */
   SCIP_CONSEXPR_EXPR*   root,                         /**< the root expression from where to start the walk */
   SCIP_DECL_CONSEXPREXPRWALK_VISIT((*enterexpr)),     /**< callback to be called when entering an expression, or NULL */
   SCIP_DECL_CONSEXPREXPRWALK_VISIT((*visitingchild)), /**< callback to be called before visiting a child, or NULL */
   SCIP_DECL_CONSEXPREXPRWALK_VISIT((*visitedchild)),  /**< callback to be called when returning from a child, or NULL */
   SCIP_DECL_CONSEXPREXPRWALK_VISIT((*leaveexpr)),     /**< callback to be called when leaving an expression, or NULL */
   void*                 data                          /**< data to be passed on to callback methods, or NULL */
   )
{
   SCIP_CONSEXPREXPRWALK_RESULT result;
   SCIP_CONSEXPR_EXPR* child;
   SCIP_CONSEXPR_EXPR* oldparent;
   int                 oldcurrentchild;

   assert(root != NULL);

   /* remember the current child of and parent of incoming root, in case we are called from within another walk */
   oldcurrentchild = root->walkcurrentchild;
   oldparent = root->walkparent;

   root->walkcurrentchild = 0;
   root->walkparent = NULL;

   /* traverse the tree */
   while( TRUE )
   {
      if( root->walkcurrentchild == 0 && enterexpr != NULL )
      {
         SCIP_CALL( (*enterexpr)(scip, root, SCIP_CONSEXPREXPRWALK_ENTEREXPR, data, &result) );
         switch( result )
         {
            case SCIP_CONSEXPREXPRWALK_CONTINUE :
               break;
            case SCIP_CONSEXPREXPRWALK_SKIP :
               root->walkcurrentchild = root->nchildren;
               break;
            case SCIP_CONSEXPREXPRWALK_ABORT :
               return SCIP_OKAY;
         }
      }

      /* if there are no more children to visit, try to go up */
      if( root->walkcurrentchild >= root->nchildren )
      {
         assert(root->walkcurrentchild == root->nchildren);

         if( leaveexpr != NULL )
         {
            SCIP_CONSEXPR_EXPR* parent = root->walkparent;   /* store parent in case the callback frees root */

            SCIP_CALL( (*leaveexpr)(scip, root, SCIP_CONSEXPREXPRWALK_LEAVEEXPR, data, &result) );
            switch( result )
            {
               case SCIP_CONSEXPREXPRWALK_CONTINUE :
                  break;
               case SCIP_CONSEXPREXPRWALK_SKIP :
                  /* this result is not allowed here */
                  SCIPABORT();
                  break;
               case SCIP_CONSEXPREXPRWALK_ABORT :
                  return SCIP_OKAY;
            }

            /* go up */
            root = parent;
         }
         else
         {
            /* go up */
            root = root->walkparent;
         }

         /* if we finished with the real root (walkparent == NULL), we are done */
         if( root == NULL )
            return SCIP_OKAY;

         /* call visitedchild of parent (now root) */
         if( visitedchild != NULL )
         {
            SCIP_CALL( (*visitedchild)(scip, root, SCIP_CONSEXPREXPRWALK_VISITEDCHILD, data, &result) );
            switch( result )
            {
               case SCIP_CONSEXPREXPRWALK_CONTINUE :
                  break;
               case SCIP_CONSEXPREXPRWALK_SKIP :
                  /* this result is not allowed here */
                  SCIPABORT();
                  break;
               case SCIP_CONSEXPREXPRWALK_ABORT :
                  return SCIP_OKAY;
            }
         }

         /* visit next child of parent (if any) */
         ++root->walkcurrentchild;

         continue;
      }

      /* keep on going down (visit next child) */
      if( visitingchild != NULL )
      {
         SCIP_CALL( (*visitingchild)(scip, root, SCIP_CONSEXPREXPRWALK_VISITINGCHILD, data, &result) );
         if( result == SCIP_CONSEXPREXPRWALK_SKIP )
         {
            /* ok, we don't go down, but skip the child: continue and try again with next child (if any) */
            ++root->walkcurrentchild;
            continue;
         }
         else if( result == SCIP_CONSEXPREXPRWALK_ABORT )
         {
            return SCIP_OKAY;
         }
      }

      /* remember the parent and set the first child that should be visited of the new root */
      child = root->children[root->walkcurrentchild];
      child->walkparent = root;
      child->walkcurrentchild = 0;
      root = child;
   }

   /* recover previous information */
   root->walkcurrentchild = oldcurrentchild;
   root->walkparent = oldparent;

   return SCIP_OKAY;
}

/** Gives the parent of an expression in an expression graph walk.
 *
 * During an expression walk, this function returns the expression from which the given expression has been accessed.
 * If not in an expression walk, the returned pointer is undefined.
 */
SCIP_CONSEXPR_EXPR* SCIPgetConsExprExprWalkParent(
   SCIP_CONSEXPR_EXPR*   expr                /**< expression which parent to get */
   )
{
   assert(expr != NULL);

   return expr->walkparent;
}

/** Gives the index of the child that will be visited next (or is currently visited) by an expression walk. */
int SCIPgetConsExprExprWalkCurrentChild(
   SCIP_CONSEXPR_EXPR*   expr                /**< expression which nextchild to get */
   )
{
   assert(expr != NULL);

   return expr->walkcurrentchild;
}

/** Gives the precedence of the expression handler of the parent expression in an expression graph walk.
 *
 * If there is no parent, then 0 is returned.
 */
unsigned int SCIPgetConsExprExprWalkParentPrecedence(
   SCIP_CONSEXPR_EXPR*   expr                /**< expression which parent to get */
   )
{
   assert(expr != NULL);

   if( expr->walkparent == NULL )
      return 0;

   return expr->walkparent->exprhdlr->precedence;
}

/*
 * constraint specific interface methods
 */

/** creates the handler for expr constraints and includes it in SCIP */
SCIP_RETCODE SCIPincludeConshdlrExpr(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONSHDLR* conshdlr;

   /* create expr constraint handler data */
   SCIP_CALL( SCIPallocClearMemory(scip, &conshdlrdata) );

   conshdlr = NULL;

   /* include constraint handler */
#if 0
   /* use SCIPincludeConshdlr() if you want to set all callbacks explicitly and realize (by getting compiler errors) when
    * new callbacks are added in future SCIP versions
    */
   SCIP_CALL( SCIPincludeConshdlr(scip, CONSHDLR_NAME, CONSHDLR_DESC,
         CONSHDLR_SEPAPRIORITY, CONSHDLR_ENFOPRIORITY, CONSHDLR_CHECKPRIORITY,
         CONSHDLR_SEPAFREQ, CONSHDLR_PROPFREQ, CONSHDLR_EAGERFREQ, CONSHDLR_MAXPREROUNDS,
         CONSHDLR_DELAYSEPA, CONSHDLR_DELAYPROP, CONSHDLR_NEEDSCONS,
         CONSHDLR_PROP_TIMING, CONSHDLR_PRESOLTIMING,
         conshdlrCopyExpr,
         consFreeExpr, consInitExpr, consExitExpr,
         consInitpreExpr, consExitpreExpr, consInitsolExpr, consExitsolExpr,
         consDeleteExpr, consTransExpr, consInitlpExpr,
         consSepalpExpr, consSepasolExpr, consEnfolpExpr, consEnfopsExpr, consCheckExpr,
         consPropExpr, consPresolExpr, consRespropExpr, consLockExpr,
         consActiveExpr, consDeactiveExpr,
         consEnableExpr, consDisableExpr, consDelvarsExpr,
         consPrintExpr, consCopyExpr, consParseExpr,
         consGetVarsExpr, consGetNVarsExpr, consGetDiveBdChgsExpr, conshdlrdata) );
#else
   /* use SCIPincludeConshdlrBasic() plus setter functions if you want to set callbacks one-by-one and your code should
    * compile independent of new callbacks being added in future SCIP versions
    */
   SCIP_CALL( SCIPincludeConshdlrBasic(scip, &conshdlr, CONSHDLR_NAME, CONSHDLR_DESC,
         CONSHDLR_ENFOPRIORITY, CONSHDLR_CHECKPRIORITY, CONSHDLR_EAGERFREQ, CONSHDLR_NEEDSCONS,
         consEnfolpExpr, consEnfopsExpr, consCheckExpr, consLockExpr,
         conshdlrdata) );
   assert(conshdlr != NULL);

   /* set non-fundamental callbacks via specific setter functions */
   SCIP_CALL( SCIPsetConshdlrActive(scip, conshdlr, consActiveExpr) );
   SCIP_CALL( SCIPsetConshdlrCopy(scip, conshdlr, conshdlrCopyExpr, consCopyExpr) );
   SCIP_CALL( SCIPsetConshdlrDeactive(scip, conshdlr, consDeactiveExpr) );
   SCIP_CALL( SCIPsetConshdlrDelete(scip, conshdlr, consDeleteExpr) );
   SCIP_CALL( SCIPsetConshdlrDelvars(scip, conshdlr, consDelvarsExpr) );
   SCIP_CALL( SCIPsetConshdlrDisable(scip, conshdlr, consDisableExpr) );
   SCIP_CALL( SCIPsetConshdlrEnable(scip, conshdlr, consEnableExpr) );
   SCIP_CALL( SCIPsetConshdlrExit(scip, conshdlr, consExitExpr) );
   SCIP_CALL( SCIPsetConshdlrExitpre(scip, conshdlr, consExitpreExpr) );
   SCIP_CALL( SCIPsetConshdlrExitsol(scip, conshdlr, consExitsolExpr) );
   SCIP_CALL( SCIPsetConshdlrFree(scip, conshdlr, consFreeExpr) );
   SCIP_CALL( SCIPsetConshdlrGetDiveBdChgs(scip, conshdlr, consGetDiveBdChgsExpr) );
   SCIP_CALL( SCIPsetConshdlrGetVars(scip, conshdlr, consGetVarsExpr) );
   SCIP_CALL( SCIPsetConshdlrGetNVars(scip, conshdlr, consGetNVarsExpr) );
   SCIP_CALL( SCIPsetConshdlrInit(scip, conshdlr, consInitExpr) );
   SCIP_CALL( SCIPsetConshdlrInitpre(scip, conshdlr, consInitpreExpr) );
   SCIP_CALL( SCIPsetConshdlrInitsol(scip, conshdlr, consInitsolExpr) );
   SCIP_CALL( SCIPsetConshdlrInitlp(scip, conshdlr, consInitlpExpr) );
   SCIP_CALL( SCIPsetConshdlrParse(scip, conshdlr, consParseExpr) );
   SCIP_CALL( SCIPsetConshdlrPresol(scip, conshdlr, consPresolExpr, CONSHDLR_MAXPREROUNDS, CONSHDLR_PRESOLTIMING) );
   SCIP_CALL( SCIPsetConshdlrPrint(scip, conshdlr, consPrintExpr) );
   SCIP_CALL( SCIPsetConshdlrProp(scip, conshdlr, consPropExpr, CONSHDLR_PROPFREQ, CONSHDLR_DELAYPROP,
         CONSHDLR_PROP_TIMING) );
   SCIP_CALL( SCIPsetConshdlrResprop(scip, conshdlr, consRespropExpr) );
   SCIP_CALL( SCIPsetConshdlrSepa(scip, conshdlr, consSepalpExpr, consSepasolExpr, CONSHDLR_SEPAFREQ, CONSHDLR_SEPAPRIORITY, CONSHDLR_DELAYSEPA) );
   SCIP_CALL( SCIPsetConshdlrTrans(scip, conshdlr, consTransExpr) );
#endif

   if( SCIPfindConshdlr(scip, "nonlinear") != NULL )
   {
      /* include the linear constraint upgrade in the linear constraint handler */
      SCIP_CALL( SCIPincludeNonlinconsUpgrade(scip, nonlinconsUpgdExpr, NULL, NONLINCONSUPGD_PRIORITY, TRUE, CONSHDLR_NAME) );
   }

   /* add expr constraint handler parameters */
   /* TODO: (optional) add constraint handler specific parameters with SCIPaddTypeParam() here */


   /* NOTE: we should not do the below when copying the constraint handler (in that case, we should call the copy callback of the expression handler */
   /* include and remember handler for variable expression */
   SCIP_CALL( SCIPincludeConsExprExprHdlrVar(scip, conshdlr) );
   assert(conshdlrdata->nexprhdlrs > 0 && strcmp(conshdlrdata->exprhdlrs[conshdlrdata->nexprhdlrs-1]->name, "var") == 0);
   conshdlrdata->exprvarhdlr = conshdlrdata->exprhdlrs[conshdlrdata->nexprhdlrs-1];

   /* include and remember handler for constant value expression */
   SCIP_CALL( SCIPincludeConsExprExprHdlrValue(scip, conshdlr) );
   assert(conshdlrdata->nexprhdlrs > 0 && strcmp(conshdlrdata->exprhdlrs[conshdlrdata->nexprhdlrs-1]->name, "val") == 0);
   conshdlrdata->exprvalhdlr = conshdlrdata->exprhdlrs[conshdlrdata->nexprhdlrs-1];

   /* include and remember handler for sum expression */
   SCIP_CALL( SCIPincludeConsExprExprHdlrSum(scip, conshdlr) );
   assert(conshdlrdata->nexprhdlrs > 0 && strcmp(conshdlrdata->exprhdlrs[conshdlrdata->nexprhdlrs-1]->name, "sum") == 0);
   conshdlrdata->exprsumhdlr = conshdlrdata->exprhdlrs[conshdlrdata->nexprhdlrs-1];

   /* include and remember handler for product expression */
   SCIP_CALL( SCIPincludeConsExprExprHdlrProduct(scip, conshdlr) );
   assert(conshdlrdata->nexprhdlrs > 0 && strcmp(conshdlrdata->exprhdlrs[conshdlrdata->nexprhdlrs-1]->name, "prod") == 0);
   conshdlrdata->exprprodhdlr = conshdlrdata->exprhdlrs[conshdlrdata->nexprhdlrs-1];

   return SCIP_OKAY;
}

/** creates and captures a expr constraint
 *
 *  @note the constraint gets captured, hence at one point you have to release it using the method SCIPreleaseCons()
 */
SCIP_RETCODE SCIPcreateConsExpr(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS**           cons,               /**< pointer to hold the created constraint */
   const char*           name,               /**< name of constraint */
   SCIP_CONSEXPR_EXPR*   expr,               /**< expression of constraint (must not be NULL) */
   SCIP_Real             lhs,                /**< left hand side of constraint */
   SCIP_Real             rhs,                /**< right hand side of constraint */
   SCIP_Bool             initial,            /**< should the LP relaxation of constraint be in the initial LP?
                                              *   Usually set to TRUE. Set to FALSE for 'lazy constraints'. */
   SCIP_Bool             separate,           /**< should the constraint be separated during LP processing?
                                              *   Usually set to TRUE. */
   SCIP_Bool             enforce,            /**< should the constraint be enforced during node processing?
                                              *   TRUE for model constraints, FALSE for additional, redundant constraints. */
   SCIP_Bool             check,              /**< should the constraint be checked for feasibility?
                                              *   TRUE for model constraints, FALSE for additional, redundant constraints. */
   SCIP_Bool             propagate,          /**< should the constraint be propagated during node processing?
                                              *   Usually set to TRUE. */
   SCIP_Bool             local,              /**< is constraint only valid locally?
                                              *   Usually set to FALSE. Has to be set to TRUE, e.g., for branching constraints. */
   SCIP_Bool             modifiable,         /**< is constraint modifiable (subject to column generation)?
                                              *   Usually set to FALSE. In column generation applications, set to TRUE if pricing
                                              *   adds coefficients to this constraint. */
   SCIP_Bool             dynamic,            /**< is constraint subject to aging?
                                              *   Usually set to FALSE. Set to TRUE for own cuts which
                                              *   are separated as constraints. */
   SCIP_Bool             removable,          /**< should the relaxation be removed from the LP due to aging or cleanup?
                                              *   Usually set to FALSE. Set to TRUE for 'lazy constraints' and 'user cuts'. */
   SCIP_Bool             stickingatnode      /**< should the constraint always be kept at the node where it was added, even
                                              *   if it may be moved to a more global node?
                                              *   Usually set to FALSE. Set to TRUE to for constraints that represent node data. */
   )
{
   /* TODO: (optional) modify the definition of the SCIPcreateConsExpr() call, if you don't need all the information */

   SCIP_CONSHDLR* conshdlr;
   SCIP_CONSDATA* consdata;

   assert(expr != NULL);

   /* find the expr constraint handler */
   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);
   if( conshdlr == NULL )
   {
      SCIPerrorMessage("expr constraint handler not found\n");
      return SCIP_PLUGINNOTFOUND;
   }

   /* create constraint data */
   SCIP_CALL( SCIPallocClearBlockMemory(scip, &consdata) );
   consdata->expr = expr;
   consdata->lhs = lhs;
   consdata->rhs = rhs;

   /* capture expression */
   SCIPcaptureConsExprExpr(consdata->expr);

   /* create constraint */
   SCIP_CALL( SCIPcreateCons(scip, cons, name, conshdlr, consdata, initial, separate, enforce, check, propagate,
         local, modifiable, dynamic, removable, stickingatnode) );

   return SCIP_OKAY;
}

/** creates and captures a expr constraint with all its constraint flags set to their
 *  default values
 *
 *  @note the constraint gets captured, hence at one point you have to release it using the method SCIPreleaseCons()
 */
SCIP_RETCODE SCIPcreateConsExprBasic(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS**           cons,               /**< pointer to hold the created constraint */
   const char*           name,               /**< name of constraint */
   SCIP_CONSEXPR_EXPR*   expr,               /**< expression of constraint (must not be NULL) */
   SCIP_Real             lhs,                /**< left hand side of constraint */
   SCIP_Real             rhs                 /**< right hand side of constraint */
   )
{
   SCIP_CALL( SCIPcreateConsExpr(scip, cons, name, expr, lhs, rhs,
         TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE) );

   return SCIP_OKAY;
}

/** returns the expression of the given expression constraint */
SCIP_CONSEXPR_EXPR* SCIPgetExprConsExpr(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint data */
   )
{
   SCIP_CONSDATA* consdata;

   assert(scip != NULL);
   assert(cons != NULL);

   if( strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(cons)), CONSHDLR_NAME) != 0 )
   {
      SCIPerrorMessage("constraint is not expression\n");
      SCIPABORT();
      return NULL;  /*lint !e527*/
   }

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   return consdata->expr;
}

/** Creates an expression from a string.
 * We specify the grammar that defines the syntax of an expression. Loosely speaking, a Base will be any "block",
 * a Factor is a Base to a power, a Term is a product of Factors and an Expression is a sum of terms
 * The actual definition:
 * <pre>
 * Expression -> ["+" | "-"] Term { ("+" | "-" | "number *") ] Term }
 * Term       -> Factor { ("*" | "/" ) Factor }
 * Factor     -> Base [ "^" "number" | "^(" "number" ")" ]
 * Base       -> "number" | "<varname>" | "(" Expression ")" | Op "(" OpExpression ")
 * </pre>
 * where [a|b] means a or b or none, (a|b) means a or b, {a} means 0 or more a.
 *
 * Note that Op and OpExpression are undefined. Op corresponds to the name of an expression handler and
 * OpExpression to whatever string the expression handler accepts (through its parse method).
 *
 * See also @ref parseExpr.
 */
SCIP_RETCODE SCIPparseConsExprExpr(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        consexprhdlr,       /**< expression constraint handler */
   const char*           exprstr,            /**< string with the expr to parse */
   const char**          finalpos,           /**< buffer to store the position of exprstr where we finished reading, or NULL if not of interest */
   SCIP_CONSEXPR_EXPR**  expr                /**< pointer to store the expr parsed */
   )
{
   const char* finalpos_;
   SCIP_RETCODE retcode;
   SCIP_HASHMAP* vartoexprvarmap;

   SCIP_CALL( SCIPhashmapCreate(&vartoexprvarmap, SCIPblkmem(scip), SCIPcalcHashtableSize(5 * SCIPgetNVars(scip))) );

   /* if parseExpr fails, we still want to free hashmap */
   retcode = parseExpr(scip, consexprhdlr, vartoexprvarmap, exprstr, &finalpos_, expr);

   SCIPhashmapFree(&vartoexprvarmap);

   if( finalpos != NULL )
      *finalpos = finalpos_;

   return retcode;
}

/** appends child to the children list of expr */
SCIP_RETCODE SCIPappendConsExprExpr(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSEXPR_EXPR*   expr,               /**< sum expression */
   SCIP_CONSEXPR_EXPR*   child               /**< expression to be appended */
   )
{
   ENSUREBLOCKMEMORYARRAYSIZE(scip, expr->children, expr->childrensize, expr->nchildren + 1);

   expr->children[expr->nchildren] = child;
   ++expr->nchildren;

   /* capture child */
   SCIPcaptureConsExprExpr(child);

   return SCIP_OKAY;
}
