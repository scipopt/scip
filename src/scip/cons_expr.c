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

/* data passed on during expression evaluation */
typedef struct
{
   SCIP_SOL*             sol;                /**< solution that is evaluated */
   unsigned int          soltag;             /**< solution tag */
   SCIP_Bool             aborted;            /**< whether the evaluation has been aborted due to an evaluation error */
} EXPREVAL_DATA;

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

/*
 * Callback methods of constraint handler
 */

/* TODO: Implement all necessary constraint handler methods. The methods with #if 0 ... #else #define ... are optional */

static
SCIP_DECL_NONLINCONSUPGD(nonlinconsUpgdExpr)
{
   assert(nupgdconss != NULL);

   *nupgdconss = 0;

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
   SCIPerrorMessage("method of expr constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

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
#if 0
static
SCIP_DECL_CONSPRINT(consPrintExpr)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of expr constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consPrintExpr NULL
#endif


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
/* Here is an attempt at defining the grammar of an expression. Will use uppercase names for variables (in the grammar sense)
 * and terminals will be between "".
 * Loosely speaking, a Base will be any "block", a Factor is a Base to a power, a Term is a product of Factors
 * and an Expression is a sum of terms
 * The actual definition:
 *
 * Expression -> ["+" | "-"] Term { ("+" | "-") Term }
 * Term       -> Factor { ("*" | "/" ) Factor }
 * Factor     -> Base [ "^" "number" | "^(" "number" ")" ]
 * Base       -> "number" | "<varname>" | "(" Expression ")" | Op "(" OpExpression ")
 *
 * where [a|b] means a or b or none, (a|b) means a or b, {a} means 0 or more a
 *
 * Note that Op and OpExpression are undefined. Op corresponds to the name of an expression handler and
 * OpExpression to whatever string the expression handler accepts (through its parse method)
 *
 * Each parse(Expr|Factor|Term|Base) returns an consexpr_expr: for the time being we assume that
 * Expr is always a exprsum
 * Term is always a exprprod
 * Factor is always a exprprod
 * Base can be anything
 *
 * @todo: we can change the grammar so that Factor becomes base and we allow a Term to be
 * Term       -> Factor { ("*" | "/" | "^") Factor }
 */
static
SCIP_RETCODE parseExpr(SCIP*, SCIP_CONSHDLR*, char* , char**, SCIP_CONSEXPR_EXPR**);

/* Base       -> "number" | "<varname>" | "(" Expression ")" | Op "(" OpExpression ")
 * builds consexprvalue, consexprvar, consexprsum or custom expr */
static
SCIP_RETCODE parseBase(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< expression constraint handler */
   char*                 expr,               /**< expr that we are parsing */
   char**                newpos,             /**< buffer to store the position of expr where we finished reading */
   SCIP_CONSEXPR_EXPR**  baseexpr            /**< buffer to store the expr parsed by Factor */
   )
{
   SCIP_VAR* var;

   //printf("parsing base from %s\n", expr);

   /* ignore whitespace */
   while( isspace((unsigned char)*expr) )
      ++expr;

   if( expr[0] == '<' )
   {
      /* parse a var */
      SCIP_CALL( SCIPparseVarName(scip, expr, &var, &expr) );

      //printf("Parsed variable %s, Should create var expression\n", SCIPvarGetName(var));
      SCIP_CALL( SCIPcreateConsExprExprVar(scip, conshdlr, baseexpr, var) );
   }
   else if( *expr == '(' )
   {
      /* parse expression */
      SCIP_CALL( parseExpr(scip, conshdlr, ++expr, newpos, baseexpr) );
      expr = *newpos;

      /* expect ')' */
      if( *expr != ')' )
      {
         SCIPerrorMessage("I was expecting a ), instead got %c from %s\n", *expr, expr);
         return SCIP_READERROR;
      }
      ++expr;
   }
   else if( isdigit(*expr) )
   {
      /* parse number */
      SCIP_Real value;
      if( !SCIPstrToRealValue(expr, &value, &expr) )
      {
         SCIPerrorMessage("error parsing number from <%s>\n", expr);
         return SCIP_READERROR;
      }
      //printf("Parsed value %g, Should create a expression value?\n", value);
      SCIP_CALL( SCIPcreateConsExprExprValue(scip, conshdlr, baseexpr, value) );
   }
   else if( isalpha(*expr) )
   {
      /* a name is coming, should find exprhandler which such name */
      char operatorname[SCIP_MAXSTRLEN];
      int i;

      /* get name */
      i = 0;
      while( *expr != '(' && !isspace((unsigned char)*expr) )
      {
         operatorname[i] = *expr;
         ++expr;
         i++;
      }
      operatorname[i] = '\0';

      /* search */
      printf("should search for exprhandler %s and give everything that is inside the parenthesis\n", operatorname);
      /* extract expression between ( ) */
      assert(expr != NULL);
   }
   else
   {
      /* Base -> "number" | "<varname>" | "(" Expression ")" | Op "(" OpExpression ") */
      SCIPerrorMessage("I was expecting a number, (exp), <varname>, Opname(Opexpr), instead got %c from %s\n", *expr, expr);
      return SCIP_READERROR;
   }
   *newpos = expr;
   return SCIP_OKAY;
}

/* Factor -> Base [ "^" "number" | "^(" "number" ")" ]
 * builds consexprprod if there is an exponent, otherwise returns base */
static
SCIP_RETCODE parseFactor(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< expression constraint handler */
   char*                 expr,               /**< expr that we are parsing */
   char**                newpos,             /**< buffer to store the position of expr where we finished reading */
   SCIP_CONSEXPR_EXPR**  factortree          /**< buffer to store the expr parsed by Factor */
   )
{
   SCIP_VAR* var;
   SCIP_CONSEXPR_EXPR*  basetree;

   //printf("parsing factor from %s\n", expr);

   /* parse Base */
   /* ignore whitespace */
   while( isspace((unsigned char)*expr) )
      ++expr;

   SCIP_CALL( parseBase(scip, conshdlr, expr, newpos, &basetree) );
   *factortree = basetree;
   expr = *newpos;

   /* check if there is a power */
   /* ignore whitespace */
   while( isspace((unsigned char)*expr) )
      ++expr;
   if( *expr == '^' )
   {
      SCIP_Real exponent;

      /* check if it is power with parenthesis */
      ++expr;
      while( isspace((unsigned char)*expr) )
         ++expr;
      if( *expr == '(' )
      {

         ++expr;
         while( isspace((unsigned char)*expr) )
            ++expr;

         /* expect '-' or a digit */
         if( isdigit(*expr) || *expr == '-' )
         {
            if( !SCIPstrToRealValue(expr, &exponent, &expr) )
            {
               SCIPerrorMessage("error parsing number from <%s>\n", expr);
               return SCIP_READERROR;
            }
         }
         else
         {
            SCIPerrorMessage("error in parsing exponent, expected '-' or a digit, received %c from <%s>\n", *expr,  expr);
            return SCIP_READERROR;
         }

         /* expect the ')' */
         while( isspace((unsigned char)*expr) )
            ++expr;
         if( *expr != ')' )
         {
            SCIPerrorMessage("error in parsing exponent: expected ')', received %c from <%s>\n", *expr,  expr);
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
            if( !SCIPstrToRealValue(expr, &exponent, &expr) )
            {
               SCIPerrorMessage("error parsing number from <%s>\n", expr);
               return SCIP_READERROR;
            }
         }
         else
         {
            SCIPerrorMessage("error in parsing exponent, expected a digit, received %c from <%s>\n", *expr,  expr);
            return SCIP_READERROR;
         }
      }

      //printf("parsed the exponen %g\n", exponent);
      /* build expression with exponent */
      SCIP_CALL( SCIPcreateConsExprExprProduct(scip, conshdlr, factortree, 1, &basetree, &exponent, 1.0) );
   }

   *newpos = expr;
   return SCIP_OKAY;
}

/* Term -> Factor { ("*" | "/" ) Factor }
 * builds consexprprod where each factor is a children
 * */
static
SCIP_RETCODE parseTerm(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< expression constraint handler */
   char*                 expr,               /**< expr that we are parsing */
   char**                newpos,             /**< buffer to store the position of expr where we finished reading */
   SCIP_CONSEXPR_EXPR**  termtree            /**< buffer to store the expr parsed by Term */
   )
{
   SCIP_CONSEXPR_EXPR* factortree;

   //printf("parsing term from %s\n", expr);

   /* parse Factor */
   /* ignore whitespace */
   while( isspace((unsigned char)*expr) )
      ++expr;

   SCIP_CALL( parseFactor(scip, conshdlr, expr, newpos, &factortree) );
   expr = *newpos;

   //printf("back to parsing Term (we have a Factor), continue parsing from %s\n", expr);
   *termtree = factortree;

   /* check if there are more factors to parse */
   while( isspace((unsigned char)*expr) )
      ++expr;
   while( *expr == '*' || *expr == '/' )
   {
      SCIP_CONSEXPR_EXPR* children[2];
      SCIP_Real exponents[2];

      exponents[0] = 1.0;
      exponents[1] = *expr == '*' ? 1.0 : -1.0;

      //printf("while parsing term, read char %c\n", *expr);
      ++expr;
      SCIP_CALL( parseFactor(scip, conshdlr, expr, newpos, &factortree) );
      expr = *newpos;

      //build (binary) tree TODO: handle power correctly, and / as well
      children[0] = *termtree;
      children[1] = factortree;
      SCIP_CALL( SCIPcreateConsExprExprProduct(scip, conshdlr, termtree, 2, children, exponents, 1.0) );

      while( isspace((unsigned char)*expr) )
         ++expr;
   }
   *newpos = expr;
   return SCIP_OKAY;
}

/* Expression -> ["+" | "-"] Term { ("+" | "-") Term }
 * builds consexprsum where each term is a children */
static
SCIP_RETCODE parseExpr(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< expression constraint handler */
   char*                 expr,               /**< expr that we are parsing */
   char**                newpos,             /**< buffer to store the position of expr where we finished reading */
   SCIP_CONSEXPR_EXPR**  exprtree            /**< buffer to store the expr parsed by Expr */
   )
{
   SCIP_Real sign;
   SCIP_CONSEXPR_EXPR* termtree;

   //printf("parsing expression %s\n", expr);

   /* ignore whitespace */
   while( isspace((unsigned char)*expr) )
      ++expr;
   // if + or -, store it
   if( *expr == '+' || *expr == '-' )
   {
      //printf("while parsing expression, read char %c\n", *expr);
      sign = *expr == '+' ? 1.0 : -1.0;
      ++expr;
   }

   SCIP_CALL( parseTerm(scip, conshdlr, expr, newpos, &termtree) );
   expr = *newpos;


   //printf("back to parsing expression (we have a term), continue parsing from %s\n", expr);
   // TODO: handle correctly the -, currently we just don't handle it
   *exprtree = termtree;
   //printf("finish parsing termtree, result:\n");
   SCIP_CALL( SCIPprintConsExprExpr(scip, termtree, NULL) );

   // loop:
   // find symbol
   // parse term
   while( isspace((unsigned char)*expr) )
      ++expr;
   while( *expr == '+' || *expr == '-' )
   {
      SCIP_CONSEXPR_EXPR* children[2]; /* just because we can't append */
      SCIP_Real coeffs[2];

      coeffs[0] = 1.0;
      coeffs[1] = *expr == '+' ? 1.0 : -1.0;

      //printf("while parsing expression, read char %c\n", *expr);
      ++expr;
      SCIP_CALL( parseTerm(scip, conshdlr, expr, newpos, &termtree) );
      expr = *newpos;

      //build (binary) tree
      children[0] = *exprtree;
      children[1] = termtree;
      //printf("building sum expr from exprtree:\n");
      //SCIP_CALL( SCIPprintConsExprExpr(scip, *exprtree, NULL) );
      //printf("and termtree:\n");
      //SCIP_CALL( SCIPprintConsExprExpr(scip, termtree, NULL) );

      SCIP_CALL( SCIPcreateConsExprExprSum(scip, conshdlr, exprtree, 2, children, coeffs, 0.0) );

      while( isspace((unsigned char)*expr) )
         ++expr;
   }
   *newpos = expr;

   return SCIP_OKAY;
}

static
SCIP_DECL_CONSPARSE(consParseExpr)
{  /*lint --e{715}*/
   SCIP_VAR** exprvars;
   SCIP_RETCODE retcode;
   int        nvars;
   SCIP_Real  lhs;
   SCIP_Real  rhs;
   const char* endptr;
   char*       nonconstendptr;
   const char* exprstart;
   const char* exprlastchar;
   int* varnames;
   int* curvarname;
   int i;
   SCIP_CONSEXPR_EXPR* consexprtree;

   SCIPdebugMessage("cons_nonlinear::consparse parsing %s\n",str);

   assert(scip != NULL);
   assert(success != NULL);
   assert(str != NULL);
   assert(name != NULL);
   assert(cons != NULL);

   /* return if string empty */
   if( !*str )
      return SCIP_OKAY;

   endptr = str;

   //expr = NULL;
   nvars = 0;

   /* set left and right hand side to their default values */
   lhs = -SCIPinfinity(scip);
   rhs =  SCIPinfinity(scip);

   /* parse constraint to get lhs, rhs, and expression in between (from cons_linear.c::consparse, but parsing whole string first, then getting expression) */

   /* check for left hand side */
   if( isdigit((unsigned char)str[0]) || ((str[0] == '-' || str[0] == '+') && isdigit((unsigned char)str[1])) )
   {
      /* there is a number coming, maybe it is a left-hand-side */
      if( !SCIPstrToRealValue(str, &lhs, &nonconstendptr) )
      {
         SCIPerrorMessage("error parsing number from <%s>\n", str);
         return SCIP_READERROR;
      }

      endptr = nonconstendptr;

      /* ignore whitespace */
      while( isspace((unsigned char)*endptr) )
         ++endptr;

      if( endptr[0] != '<' || endptr[1] != '=' )
      {
         /* no '<=' coming, so it was the first coefficient, but not a left-hand-side */
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

   //printf("str should start at beginning of expr: %s\n", str);

   /* Move endptr forward until we find end of expression */
   endptr = str;
   while( !(strncmp(endptr, "[free]", 6) == 0)    &&
          !(endptr[0] == '<' && endptr[1] == '=') &&
          !(endptr[0] == '=' && endptr[1] == '=') &&
          !(endptr[0] == '>' && endptr[1] == '=') &&
          !(endptr[0] == '\0') )
      ++endptr;

   //printf("just found end of expr: %s\n", endptr);
   exprstart = str;
   exprlastchar = endptr - 1;

   *success = FALSE;
   str = endptr;

   /* check for left or right hand side */
   while( isspace((unsigned char)*str) )
      ++str;

   /* check for free constraint */
   if( strncmp(str, "[free]", 6) == 0 )
   {
      if( !SCIPisInfinity(scip, -lhs) )
      {
         SCIPerrorMessage("cannot have left hand side and [free] status \n");
         return SCIP_OKAY;
      }
      (*success) = TRUE;
   }
   else
   {
      switch( *str )
      {
         case '<':
            *success = SCIPstrToRealValue(str+2, &rhs, &nonconstendptr);
            break;
         case '=':
            if( !SCIPisInfinity(scip, -lhs) )
            {
               SCIPerrorMessage("cannot have == on rhs if there was a <= on lhs\n");
               return SCIP_OKAY;
            }
            else
            {
               *success = SCIPstrToRealValue(str+2, &rhs, &nonconstendptr);
               lhs = rhs;
            }
            break;
         case '>':
            if( !SCIPisInfinity(scip, -lhs) )
            {
               SCIPerrorMessage("cannot have => on rhs if there was a <= on lhs\n");
               return SCIP_OKAY;
            }
            else
            {
               *success = SCIPstrToRealValue(str+2, &lhs, &nonconstendptr);
               break;
            }
         case '\0':
            *success = TRUE;
            break;
         default:
            SCIPerrorMessage("unexpected character %c\n", *str);
            return SCIP_OKAY;
      }
   }

   //printf("This are the current string before start parsing:\nstr: %s\nexprstart: %s\nexprlastchar: %s\nendptr: %s\n", str, exprstart, exprlastchar, endptr);

   /* alloc some space for variable names incl. indices; shouldn't be longer than expression string, and we even give it sizeof(int) times this length (plus 5) */
   SCIP_CALL( SCIPallocBufferArray(scip, &varnames, (int) (exprlastchar - exprstart) + 5) );

   /* parse expression */
   {
      char newstr[SCIP_MAXSTRLEN];
      char* newpos;

      snprintf(newstr, exprlastchar - exprstart + 1, "%s", exprstart);
      printf("here i should start parsing %s\n", newstr);
      /* should probably process the return value separately to free data when fail */
      SCIP_CALL( parseExpr(scip, conshdlr, newstr, &newpos, &consexprtree) );

      printf("finish parsing, printing what i got:\n");
      SCIPinfoMessage(scip, NULL, "\n");
      SCIP_CALL( SCIPprintConsExprExpr(scip, consexprtree, NULL) );
      SCIPinfoMessage(scip, NULL, "\n");
      printf("done printing\n");
   }
   //retcode = SCIPexprParse(SCIPblkmem(scip), SCIPgetMessagehdlr(scip), &expr, exprstart, exprlastchar, &nvars, varnames);

   /* create constraint */
   SCIP_CALL( SCIPcreateConsExpr(scip, cons, name,
      consexprtree, lhs, rhs,
      initial, separate, enforce, check, propagate, local, modifiable, dynamic, removable, stickingatnode) );
   assert(cons != NULL);
   assert(*cons != NULL);

   //printf("created expression constraint: <%s>\n", SCIPconsGetName(*cons));
   //SCIP_CALL( SCIPprintConsExprExpr(scip, (*cons)->expr, NULL) );
   //SCIPinfoMessage(scip, NULL, "\n");

   //SCIP_CALL( SCIPexprtreeFree(&exprtree) );

 TERMINATE:
   //SCIPfreeBufferArrayNull(scip, &exprvars);
   //SCIPfreeBufferArrayNull(scip, &varnames);

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
   int                         precedence,   /**< precedence of expression operation (used for printing) */
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
int SCIPgetConsExprExprHdlrPrecedence(
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
   SCIPinfoMessage(scip, NULL, "\n");

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

/** gives the value from the last evaluation of an expression (or SCIP_INVALID if there was an eval error) */
SCIP_Real SCIPgetConsExprExprValue(
   SCIP_CONSEXPR_EXPR*     expr              /**< expression */
   )
{
   assert(expr != NULL);

   return expr->evalvalue;
}

/** gives the evaluation tag from the last evaluation, or 0 */
unsigned int SCIPgetConsExprExprEvalTag(
   SCIP_CONSEXPR_EXPR*     expr              /**< expression */
   )
{
   assert(expr != NULL);

   return expr->evaltag;
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

/** walks the expression graph in depth-first manner and executes callbacks at certain places
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

   assert(root != NULL);

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
   return expr->walkparent;
}

/** Gives the index of the child that will be visited next (or is currently visited) by an expression walk. */
int SCIPgetConsExprExprWalkCurrentChild(
   SCIP_CONSEXPR_EXPR*   expr                /**< expression which nextchild to get */
   )
{
   return expr->walkcurrentchild;
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
