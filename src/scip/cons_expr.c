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
   SCIP_CONSEXPR_EXPR*         expr;        /**< expression that represents this constraint (must evaluate to 0 (FALSE) or 1 (TRUE)) */
};

/** constraint handler data */
struct SCIP_ConshdlrData
{
   SCIP_CONSEXPR_OPERANDHDLR** ophdlrs;     /**< operand handlers */
   int                         nophdlrs;    /**< number of operand handlers */
   int                         ophdlrssize; /**< size of ophdlrs array */

   SCIP_CONSEXPR_OPERANDHDLR*  opvarhdlr;    /**< variable operand handler */
   SCIP_CONSEXPR_OPERANDHDLR*  opvalhdlr;    /**< value operand handler */
   SCIP_CONSEXPR_OPERANDHDLR*  opsumhdlr;    /**< summation operand handler */
   SCIP_CONSEXPR_OPERANDHDLR*  opprodhdlr;   /**< product operand handler */
};


/*
 * Local methods
 */

/* put your local methods here, and declare them static */

static
SCIP_RETCODE freeExpr(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSEXPR_EXPR**  expr                /**< expression to be freed */
   )
{
   SCIP_CONSEXPR_EXPR** children;
   int c = 0;

   assert(expr != NULL);
   assert(*expr != NULL);

   /* free operator data */
   SCIP_CALL( SCIPfreeOperandData(scip, (*expr)->ophdlr, &(*expr)->opdata, SCIPgetConsExprExprNChildren(*expr)) );

   /* release children */
   children = SCIPgetConsExprExprChildren(*expr);
   for( c = SCIPgetConsExprExprNChildren(*expr)-1; c >= 0; --c )
   {
      SCIP_CALL( SCIPreleaseConsExprExpr(scip, &children[c]) );
   }

   /* free children array, if multivariate */
   if( (*expr)->variability == SCIP_CONSEXPR_MULTIVARIATE )
   {
      SCIPfreeBlockMemoryArrayNull(scip, (*expr)->children.array.children, (*expr)->children.array.childrensize);
   }

   SCIPfreeBlockMemory(scip, expr);
   assert(*expr == NULL);

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
   SCIP_CONSEXPR_OPERANDHDLR* ophdlr;
   int i;

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   for( i = 0; i < conshdlrdata->nophdlrs; ++i )
   {
      ophdlr = conshdlrdata->ophdlrs[i];
      assert(ophdlr != NULL);

      if( ophdlr->freehdlr != NULL )
      {
         SCIP_CALL( (*ophdlr->freehdlr)(scip, conshdlr, ophdlr, &ophdlr->data) );
      }

      SCIPfreeMemory(scip, &ophdlr->name);
      SCIPfreeMemoryNull(scip, &ophdlr->desc);

      SCIPfreeMemory(scip, &ophdlr);
   }

   SCIPfreeBlockMemoryArray(scip, &conshdlrdata->ophdlrs, conshdlrdata->ophdlrssize);

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
#if 0
static
SCIP_DECL_CONSDELETE(consDeleteExpr)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of expr constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consDeleteExpr NULL
#endif


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
#if 0
static
SCIP_DECL_CONSPARSE(consParseExpr)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of expr constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consParseExpr NULL
#endif


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



/** creates the handler for an expression operand and includes it into the expression constraint handler */
SCIP_RETCODE SCIPincludeOperandHdlrBasic(
   SCIP*                      scip,          /**< SCIP data structure */
   SCIP_CONSHDLR*             conshdlr,      /**< expression constraint handler */
   SCIP_CONSEXPR_OPERANDHDLR** ophdlr,       /**< buffer where to store operand handler */
   const char*                name,          /**< name of operand (must not be NULL) */
   const char*                desc,          /**< description of operand (can be NULL) */
   SCIP_CONSEXPR_OPERANDHDLRDATA* data       /**< data of operand handler */
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata;

   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);
   assert(name != NULL);
   assert(ophdlr != NULL);

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   SCIP_CALL( SCIPallocClearMemory(scip, ophdlr) );

   SCIP_CALL( SCIPduplicateMemoryArray(scip, &(*ophdlr)->name, name, strlen(name)+1) );
   if( desc != NULL )
   {
      SCIP_CALL( SCIPduplicateMemoryArray(scip, &(*ophdlr)->desc, desc, strlen(desc)+1) );
   }

   (*ophdlr)->data = data;

   ENSUREBLOCKMEMORYARRAYSIZE(scip, conshdlrdata->ophdlrs, conshdlrdata->ophdlrssize, conshdlrdata->nophdlrs+1);

   conshdlrdata->ophdlrs[conshdlrdata->nophdlrs] = *ophdlr;
   ++conshdlrdata->nophdlrs;

   return SCIP_OKAY;
}

/** set the operand handler callbacks to copy and free an operand handler */
SCIP_RETCODE SCIPsetOperandHdlrCopyFreeHdlr(
   SCIP*                      scip,          /**< SCIP data structure */
   SCIP_CONSHDLR*             conshdlr,      /**< expression constraint handler */
   SCIP_CONSEXPR_OPERANDHDLR* ophdlr,        /**< buffer where to store operand handler */
   SCIP_DECL_CONSEXPR_OPERANDCOPYHDLR((*copyhdlr)), /**< handler copy method (can be NULL) */
   SCIP_DECL_CONSEXPR_OPERANDFREEHDLR((*freehdlr)) /**< handler free method (can be NULL) */
)
{
   assert(ophdlr != NULL);

   ophdlr->copyhdlr = copyhdlr;
   ophdlr->freehdlr = freehdlr;

   return SCIP_OKAY;
}

/** set the operand handler callbacks to copy and free operand data */
SCIP_RETCODE SCIPsetOperandHdlrCopyFreeData(
   SCIP*                      scip,          /**< SCIP data structure */
   SCIP_CONSHDLR*             conshdlr,      /**< expression constraint handler */
   SCIP_CONSEXPR_OPERANDHDLR* ophdlr,        /**< buffer where to store operand handler */
   SCIP_DECL_CONSEXPR_OPERANDCOPYDATA((*copydata)), /**< copy method of operand data (can be NULL for operands without data) */
   SCIP_DECL_CONSEXPR_OPERANDFREEDATA((*freedata))  /**< free method of operand data (can be NULL if data does not need to be freed) */
)
{
   assert(ophdlr != NULL);

   ophdlr->copydata = copydata;
   ophdlr->freedata = freedata;

   return SCIP_OKAY;
}

/** set the print callback of an operand handler */
SCIP_RETCODE SCIPsetOperandHdlrPrint(
   SCIP*                      scip,          /**< SCIP data structure */
   SCIP_CONSHDLR*             conshdlr,      /**< expression constraint handler */
   SCIP_CONSEXPR_OPERANDHDLR* ophdlr,        /**< buffer where to store operand handler */
   SCIP_DECL_CONSEXPR_OPERANDPRINT((*print)) /**< print method of operand data (can be NULL) */
)
{
   assert(ophdlr != NULL);

   ophdlr->print = print;

   return SCIP_OKAY;
}

/** gives the name of an operand handler */
const char* SCIPgetOperandHdlrName(
   SCIP_CONSEXPR_OPERANDHDLR* ophdlr         /**< operand handler */
)
{
   assert(ophdlr != NULL);

   return ophdlr->name;
}

/** gives the description of an operand handler (can be NULL) */
const char* SCIPgetOperandHdlrDescription(
   SCIP_CONSEXPR_OPERANDHDLR* ophdlr         /**< operand handler */
)
{
   assert(ophdlr != NULL);

   return ophdlr->desc;
}

/** gives the data of an operand handler */
SCIP_CONSEXPR_OPERANDHDLRDATA* SCIPgetOperandHdlrData(
   SCIP_CONSEXPR_OPERANDHDLR* ophdlr         /**< operand handler */
)
{
   assert(ophdlr != NULL);

   return ophdlr->data;
}

/** frees operand data */
SCIP_RETCODE SCIPfreeOperandData(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSEXPR_OPERANDHDLR* ophdlr,        /**< operand handler */
   SCIP_CONSEXPR_OPERANDDATA** opdata,       /**< operand data to be freed, if *opdata is not NULL */
   int                   nchildren           /**< number of children of corresponding expression */
)
{
   assert(opdata != NULL);

   /* if nothing to free, then free nothing */
   if( *opdata == NULL )
      return SCIP_OKAY;

   /* if there is no callback to free data, then free nothing */
   if( ophdlr->freedata == NULL )
   {
      *opdata = NULL;
      return SCIP_OKAY;
   }

   SCIP_CALL( ophdlr->freedata(scip, ophdlr, opdata, nchildren) );
   assert(*opdata == NULL);

   return SCIP_OKAY;
}



/** creates and captures a multivariate expression with given operand and children */
SCIP_RETCODE SCIPcreateConsExprExprMultivariate(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        consexprhdlr,       /**< expression constraint handler */
   SCIP_CONSEXPR_EXPR**  expr,               /**< pointer where to store expression */
   SCIP_CONSEXPR_OPERANDHDLR* ophdlr,        /**< operand handler */
   SCIP_CONSEXPR_OPERANDDATA* opdata,        /**< operand data */
   int                   nchildren,          /**< number of children */
   SCIP_CONSEXPR_EXPR*   children            /**< children */
   )
{
   int c;

   assert(expr != NULL);
   assert(ophdlr != NULL);
   assert(children != NULL || nchildren == 0);

   SCIP_CALL( SCIPallocBlockMemory(scip, expr) );

   (*expr)->ophdlr = ophdlr;
   (*expr)->opdata = opdata;
   (*expr)->variability = SCIP_CONSEXPR_MULTIVARIATE;

   if( nchildren > 0 )
   {
      SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &(*expr)->children.array.children, children, nchildren) );
      (*expr)->children.array.nchildren = nchildren;
      (*expr)->children.array.childrensize = nchildren;

      for( c = 0; c < nchildren; ++c )
         SCIPcaptureConsExprExpr((*expr)->children.array.children[c]);
   }
   else
   {
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(*expr)->children.array.children, 2) );
      (*expr)->children.array.nchildren = 0;
      (*expr)->children.array.childrensize = 2;
   }

   (*expr)->nuses = 0;
   SCIPcaptureConsExprExpr(*expr);

   return SCIP_OKAY;
}

/** creates and captures a bivariate expression with given operand and children */
SCIP_RETCODE SCIPcreateConsExprExprBivariate(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        consexprhdlr,       /**< expression constraint handler */
   SCIP_CONSEXPR_EXPR**  expr,               /**< pointer where to store expression */
   SCIP_CONSEXPR_OPERANDHDLR* ophdlr,        /**< operand handler */
   SCIP_CONSEXPR_OPERANDDATA* opdata,        /**< operand data */
   SCIP_CONSEXPR_EXPR*   child1,             /**< first child */
   SCIP_CONSEXPR_EXPR*   child2              /**< second child */
   )
{
   assert(expr != NULL);
   assert(ophdlr != NULL);
   assert(child1 != NULL);
   assert(child2 != NULL);

   SCIP_CALL( SCIPallocBlockMemory(scip, expr) );

   (*expr)->ophdlr = ophdlr;
   (*expr)->opdata = opdata;
   (*expr)->variability = SCIP_CONSEXPR_BIVARIATE;
   (*expr)->children.pair[0] = child1;
   (*expr)->children.pair[1] = child2;

   SCIPcaptureConsExprExpr(child1);
   SCIPcaptureConsExprExpr(child2);

   (*expr)->nuses = 0;
   SCIPcaptureConsExprExpr(*expr);

   return SCIP_OKAY;
}

/** creates and captures a univariate expression with given operand and child */
SCIP_RETCODE SCIPcreateConsExprExprUnivariate(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        consexprhdlr,       /**< expression constraint handler */
   SCIP_CONSEXPR_EXPR**  expr,               /**< pointer where to store expression */
   SCIP_CONSEXPR_OPERANDHDLR* ophdlr,        /**< operand handler */
   SCIP_CONSEXPR_OPERANDDATA* opdata,        /**< operand data */
   SCIP_CONSEXPR_EXPR*   child               /**< child */
   )
{
   assert(expr != NULL);
   assert(ophdlr != NULL);
   assert(child != NULL);

   SCIP_CALL( SCIPallocBlockMemory(scip, expr) );

   (*expr)->ophdlr = ophdlr;
   (*expr)->opdata = opdata;
   (*expr)->variability = SCIP_CONSEXPR_UNIVARIATE;
   (*expr)->children.single = child;

   SCIPcaptureConsExprExpr(child);

   (*expr)->nuses = 0;
   SCIPcaptureConsExprExpr(*expr);

   return SCIP_OKAY;
}

/** creates and captures a variate expression with given operand */
SCIP_RETCODE SCIPcreateConsExprExprInvariate(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        consexprhdlr,       /**< expression constraint handler */
   SCIP_CONSEXPR_EXPR**  expr,               /**< pointer where to store expression */
   SCIP_CONSEXPR_OPERANDHDLR* ophdlr,        /**< operand handler */
   SCIP_CONSEXPR_OPERANDDATA* opdata         /**< operand data */
   )
{
   assert(expr != NULL);
   assert(ophdlr != NULL);

   SCIP_CALL( SCIPallocBlockMemory(scip, expr) );

   (*expr)->ophdlr = ophdlr;
   (*expr)->opdata = opdata;
   (*expr)->variability = SCIP_CONSEXPR_INVARIATE;

   (*expr)->nuses = 0;
   SCIPcaptureConsExprExpr(*expr);

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
      SCIP_CALL( freeExpr(scip, expr) );

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

   if( expr->variability == SCIP_CONSEXPR_MULTIVARIATE )
      return expr->children.array.nchildren;

   /* the enum values correspond to the number of children, if not multivariate */
   return (int)expr->variability;
}

/** gives the child of a univariate expression */
SCIP_CONSEXPR_EXPR* SCIPgetConsExprExprChild(
   SCIP_CONSEXPR_EXPR*   expr               /**< expression */
   )
{
   assert(expr != NULL);
   assert(expr->variability == SCIP_CONSEXPR_UNIVARIATE);

   return expr->children.single;
}

/** gives the children of a non-invariate expression */
SCIP_CONSEXPR_EXPR** SCIPgetConsExprExprChildren(
   SCIP_CONSEXPR_EXPR*   expr               /**< expression */
   )
{
   assert(expr != NULL);
   assert(expr->variability != SCIP_CONSEXPR_INVARIATE);

   switch( expr->variability )
   {
      case SCIP_CONSEXPR_UNIVARIATE :
         return &expr->children.single;
      case SCIP_CONSEXPR_BIVARIATE :
         return expr->children.pair;
      case SCIP_CONSEXPR_MULTIVARIATE :
         return expr->children.array.children;
      default:
         /* invariate case handled in assert above */
         SCIPABORT();
         return NULL;
   }
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


   /* NOTE: we should not do the below when copying the constraint handler (in that case, we should call the copy callback of the operator handler */
   /* include and remember handler for variable operator */
   SCIP_CALL( SCIPincludeOperandHdlrVar(scip, conshdlr) );
   assert(conshdlrdata->nophdlrs > 0 && strcmp(conshdlrdata->ophdlrs[conshdlrdata->nophdlrs-1]->name, "var") == 0);
   conshdlrdata->opvarhdlr = conshdlrdata->ophdlrs[conshdlrdata->nophdlrs-1];

   /* include and remember handler for constant value operator */
   SCIP_CALL( SCIPincludeOperandHdlrValue(scip, conshdlr) );
   assert(conshdlrdata->nophdlrs > 0 && strcmp(conshdlrdata->ophdlrs[conshdlrdata->nophdlrs-1]->name, "val") == 0);
   conshdlrdata->opvalhdlr = conshdlrdata->ophdlrs[conshdlrdata->nophdlrs-1];

   /* include and remember handler for sum operator */
   SCIP_CALL( SCIPincludeOperandHdlrSum(scip, conshdlr) );
   assert(conshdlrdata->nophdlrs > 0 && strcmp(conshdlrdata->ophdlrs[conshdlrdata->nophdlrs-1]->name, "sum") == 0);
   conshdlrdata->opsumhdlr = conshdlrdata->ophdlrs[conshdlrdata->nophdlrs-1];

   /* include and remember handler for product operator */
   SCIP_CALL( SCIPincludeOperandHdlrProduct(scip, conshdlr) );
   assert(conshdlrdata->nophdlrs > 0 && strcmp(conshdlrdata->ophdlrs[conshdlrdata->nophdlrs-1]->name, "prod") == 0);
   conshdlrdata->opsumhdlr = conshdlrdata->ophdlrs[conshdlrdata->nophdlrs-1];

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
   int                   nvars,              /**< number of variables in the constraint */
   SCIP_VAR**            vars,               /**< array with variables of constraint entries */
   SCIP_Real*            coefs,              /**< array with coefficients of constraint entries */
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

   SCIPerrorMessage("method of expr constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527} --e{715}*/

   /* find the expr constraint handler */
   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);
   if( conshdlr == NULL )
   {
      SCIPerrorMessage("expr constraint handler not found\n");
      return SCIP_PLUGINNOTFOUND;
   }

   /* create constraint data */
   consdata = NULL;
   /* TODO: create and store constraint specific data here */

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
   int                   nvars,              /**< number of variables in the constraint */
   SCIP_VAR**            vars,               /**< array with variables of constraint entries */
   SCIP_Real*            coefs,              /**< array with coefficients of constraint entries */
   SCIP_Real             lhs,                /**< left hand side of constraint */
   SCIP_Real             rhs                 /**< right hand side of constraint */
   )
{
   SCIP_CALL( SCIPcreateConsExpr(scip, cons, name, nvars, vars, coefs, lhs, rhs,
         TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE) );

   return SCIP_OKAY;
}
