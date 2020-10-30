/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2020 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not visit scip.zib.de.         */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   expr.c
 * @brief  functions for algebraic expressions
 * @author Ksenia Bestuzheva
 * @author Benjamin Mueller
 * @author Felipe Serrano
 * @author Stefan Vigerske
 */

/* #define PARSE_DEBUG */

#include <string.h>
#include <ctype.h>
#include <assert.h>

#include "scip/scip_expr.h"
#include "scip/expr.h"
#include "scip/struct_expr.h"
#include "nlpi/nlpi_ipopt.h" /* for LAPACK */

/*
 * Data structures
 */

/** variable mapping data passed on during copying expressions when copying SCIP instances */
typedef struct
{
   SCIP_HASHMAP*         varmap;             /**< SCIP_HASHMAP mapping variables of the source SCIP to corresponding variables of the target SCIP */
   SCIP_HASHMAP*         consmap;            /**< SCIP_HASHMAP mapping constraints of the source SCIP to corresponding constraints of the target SCIP */
   SCIP_Bool             global;             /**< should a global or a local copy be created */
   SCIP_Bool             valid;              /**< indicates whether every variable copy was valid */
} COPY_MAPVAR_DATA;

/** printing to file data */
struct SCIP_ExprPrintData
{
   FILE*                 file;               /**< file to print to */
   SCIP_EXPRITER*        iterator;           /**< iterator to use */
   SCIP_Bool             closefile;          /**< whether file need to be closed when finished printing */
   SCIP_HASHMAP*         leaveexprs;         /**< hashmap storing leave (no children) expressions */
   SCIP_EXPRPRINT_WHAT   whattoprint;        /**< flags that indicate what to print for each expression */
};

/*
 * Local methods
 */

/** creates an expression */
static
SCIP_RETCODE createExpr(
   SCIP_SET*             set,                /**< global SCIP settings */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_EXPR**           expr,               /**< pointer where to store expression */
   SCIP_EXPRHDLR*        exprhdlr,           /**< expression handler */
   SCIP_EXPRDATA*        exprdata,           /**< expression data (expression assumes ownership) */
   int                   nchildren,          /**< number of children */
   SCIP_EXPR**           children,           /**< children (can be NULL if nchildren is 0) */
   )
{
   int c;

   assert(set != NULL);
   assert(blkmem != NULL);
   assert(expr != NULL);
   assert(exprhdlr != NULL);
   assert(children != NULL || nchildren == 0);
   assert(exprdata == NULL || exprhdlr->copydata != NULL); /* copydata must be available if there is expression data */
   assert(exprdata == NULL || exprhdlr->freedata != NULL); /* freedata must be available if there is expression data */

   SCIP_ALLOC( BMSallocClearBlockMemory(blkmem, expr) );

   (*expr)->exprhdlr = exprhdlr;
   (*expr)->exprdata = exprdata;
   (*expr)->curvature = SCIP_EXPRCURV_UNKNOWN;

   /* initialize activity to entire interval */
   SCIPintervalSetEntire(SCIP_INTERVAL_INFINITY, &(*expr)->activity);

   if( nchildren > 0 )
   {
      SCIP_ALLOC( BMSduplicateBlockMemoryArray(blkmem, &(*expr)->children, children, nchildren) );
      (*expr)->nchildren = nchildren;
      (*expr)->childrensize = nchildren;

      for( c = 0; c < nchildren; ++c )
         SCIPcaptureExpr((*expr)->children[c]);
   }

   SCIPcaptureExpr(*expr);

   return SCIP_OKAY;
}

/** initializes the ownerdata of an expression
 *
 * typically called right after creating an expression
 */
static
SCIP_RETCODE createExprOwnerData(
   SCIP_SET*             set,                          /**< global SCIP settings */
   SCIP_EXPR*            expr,                         /**< expression for which to create ownerdata */
   SCIP_DECL_EXPR_OWNERDATACREATE((*ownerdatacreate)), /**< function to call to create ownerdata */
   SCIP_EXPR_OWNERDATACREATEDATA* ownerdatacreatedata, /**< data to pass to ownerdatacreate */
   SCIP_DECL_EXPR_OWNERDATAFREE((*ownerdatafree))      /**< function to call when freeing expression, e.g., to free ownerdata */
   )
{
   assert(set != NULL);
   assert(expr != NULL);

   /* expr should not yet have ownerdata or ownerdatafree
    * (if this becomes an issue some day, we could call the ownerdatafree here instead of the asserts)
    */
   assert(expr->ownerdata == NULL);
   assert(expr->ownerdatafree == NULL);

   if( ownerdatacreate != NULL )
   {
      SCIP_CALL( ownerdatacreate(set->scip, expr, &expr->ownerdata, ownerdatacreatedata) );
   }
   expr->ownerdatafree = ownerdatafree;

   return SCIP_OKAY;
}

/** frees an expression */
static
SCIP_RETCODE freeExpr(
   SCIP_SET*             set,                /**< global SCIP settings */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_EXPR**           expr                /**< pointer to free the expression */
   )
{
   assert(scip != NULL);
   assert(expr != NULL);
   assert(*expr != NULL);
   assert((*expr)->nuses == 1);
   assert((*expr)->quaddata == NULL);

   /* call ownerdatafree callback, if given
    * we intentially call this also if ownerdata is NULL, so owner can be notified without storing data
    */
   if( (*expr)->ownerdatafree != NULL )
   {
      SCIP_CALL( (*expr)->ownerdatafree(set->scip, *expr, &(*expr)->ownerdata) );
   }
   assert((*expr)->ownerdata == NULL);

   /* free children array, if any */
   BMSfreeBlockMemoryArrayNull(blkmem, &(*expr)->children, (*expr)->childrensize);

   BMSfreeBlockMemory(blkmem, expr);
   assert(*expr == NULL);

   return SCIP_OKAY;
}

/** variable mapping callback to call when copying expressions (within same or different SCIPs) */
static
SCIP_DECL_EXPR_MAPVAR(copyVar)
{
   COPY_MAPVAR_DATA* data;
   SCIP_Bool valid;

   assert(sourcevar != NULL);
   assert(targetvar != NULL);
   assert(mapvardata != NULL);

   data = (COPY_MAPVAR_DATA*)mapvardata;

   SCIP_CALL( SCIPgetVarCopy(sourcescip, targetscip, sourcevar, targetvar, data->varmap, data->consmap, data->global, &valid) );
   assert(*targetvar != NULL);

   /* if copy was not valid, store so in mapvar data */
   if( !valid )
      data->valid = FALSE;

   /* caller assumes that target variable has been captured */
   SCIP_CALL( SCIPcaptureVar(targetscip, *targetvar) );

   return SCIP_OKAY;
}

/** copies an expression including subexpressions
 *
 * @note If copying fails due to an expression handler not being available in the targetscip, then *targetexpr will be set to NULL.
 *
 * Variables can be mapped to different ones by specifying a mapvar callback.
 * For all or some expressions, a mapping to an existing expression can be specified via the mapexpr callback.
 * The mapped expression (including its children) will not be copied in this case and its ownerdata will not be touched.
 * If, however, the mapexpr callback returns NULL for the targetexpr, then the expr will be copied in the usual way.
 */
static
SCIP_RETCODE copyExpr(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< dynamic problem statistics */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             targetset,          /**< global SCIP settings data structure where target expression will live */
   BMS_BLKMEM*           targetblkmem,       /**< block memory in target SCIP */
   SCIP_EXPR*            sourceexpr,         /**< expression to be copied */
   SCIP_EXPR**           targetexpr,         /**< buffer to store pointer to copy of source expression */
   SCIP_DECL_EXPR_MAPVAR((*mapvar)),         /**< variable mapping function, or NULL for identity mapping */
   void*                 mapvardata,         /**< data of variable mapping function */
   SCIP_DECL_EXPR_MAPEXPR((*mapexpr)),       /**< expression mapping function, or NULL for creating new expressions */
   void*                 mapexprdata,        /**< data of expression mapping function */
   SCIP_DECL_EXPR_OWNERDATACREATE((*ownerdatacreate)), /**< function to call on expression copy to create ownerdata */
   SCIP_EXPR_OWNERDATACREATEDATA* ownerdatacreatedata, /**< data to pass to ownerdatacreate */
   SCIP_DECL_EXPR_OWNERDATAFREE((*ownerdatafree)),     /**< function to call when freeing expression, e.g., to free ownerdata */
   )
{
   SCIP_EXPRITER* it;
   SCIP_EXPRITER_USERDATA expriteruserdata;
   SCIP_EXPR* expr;
   SCIP* sourcescip = set->scip;        /* SCIP data structure corresponding to source expression */
   SCIP* targetscip = targetset->scip;  /* SCIP data structure where target expression will live */

   assert(set != NULL);
   assert(stat != NULL);
   assert(blkmem != NULL);
   assert(targetset != NULL);
   assert(sourceexpr != NULL);
   assert(targetexpr != NULL);
   assert(sourcescip != NULL);
   assert(targetscip != NULL);

   SCIP_CALL( SCIPexpriteratorCreate(stat, blkmem, &it) );
   SCIP_CALL( SCIPexpriteratorInit(it, sourceexpr, SCIP_EXPRITER_DFS, TRUE) );  /*TODO use FALSE, i.e., don't duplicate common subexpr? */
   SCIPexpriteratorSetStagesDFS(it, SCIP_EXPRITER_ENTEREXPR | SCIP_EXPRITER_VISITEDCHILD);

   expr = sourceexpr;
   while( !SCIPexpriteratorIsEnd(it) )
   {
      switch( SCIPexpriteratorGetStageDFS(it) )
      {
         case SCIP_EXPRITER_ENTEREXPR :
         {
            /* create expr that will hold the copy */
            SCIP_EXPR* exprcopy = NULL;

            if( mapexprdata != NULL )
            {
               SCIP_CALL( mapexprdata(targetscip, &exprcopy, sourcescip, expr, mapexprdata) );
               if( exprcopy != NULL )
               {
                  /* map callback gave us an expression to use for the copy */
                  /* store targetexpr */
                  expriteruserdata.ptrval = exprcopy;
                  SCIPexpriteratorSetCurrentUserData(it, expriteruserdata);

                  /* skip subexpression (assume that exprcopy is a complete copy) and continue */
                  expr = SCIPexpriteratorSkipDFS(it);
                  continue;
               }
            }

            /* if the source is a variable expression create a variable expression directly; otherwise copy the expression data */
            if( SCIPisExprVar(expr) )
            {
               SCIP_VAR* sourcevar;
               SCIP_VAR* targetvar;

               sourcevar = SCIPgetConsExprExprVarVar(expr);
               assert(sourcevar != NULL);
               targetvar = NULL;

               /* get the corresponding variable in the target SCIP */
               if( mapvar != NULL )
               {
                  SCIP_CALL( mapvar(targetscip, &targetvar, sourcescip, sourcevar, mapvardata) );
                  SCIP_CALL( SCIPcreateConsExprExprVar(targetscip, &exprcopy, targetvar) );

                  /* we need to release once since it has been captured by the mapvar() and createExprVar() call */
                  SCIP_CALL( SCIPreleaseVar(targetscip, &targetvar) );
               }
               else
               {
                  targetvar = sourcevar;
                  SCIP_CALL( SCIPcreateConsExprExprVar(targetscip, &exprcopy, targetvar) );
               }
            }
            else
            {
               SCIP_EXPRHDLR* targetexprhdlr;
               SCIP_EXPRDATA* targetexprdata;

               /* get the exprhdlr of the target scip */
               if( targetscip != sourcescip )
               {
                  targetexprhdlr = SCIPsetFindExprhdlr(targetset, SCIPexprhdlrGetName(SCIPexprGetHdlr(expr)));

                  if( targetexprhdlr == NULL )
                  {
                     /* expression handler not in target scip (probably did not have a copy callback) -> abort */
                     expriteruserdata.ptrval = NULL;
                     SCIPexpriteratorSetCurrentUserData(it, expriteruserdata);

                     expr = SCIPexpriteratorSkipDFS(it);
                     continue;
                  }
               }
               else
               {
                  targetexprhdlr = SCIPexprGetHdlr(expr);
               }
               assert(targetexprhdlr != NULL);

               /* copy expression data */
               if( expr->exprdata != NULL )
               {
                  assert(expr->exprhdlr->copydata != NULL);
                  SCIP_CALL( expr->exprhdlr->copydata(targetscip, targetexprhdlr, &targetexprdata, sourcescip, expr, mapvar, mapvardata) );
               }
               else
               {
                  targetexprdata = NULL;
               }

               /* create in targetexpr an expression of the same type as expr, but without children for now */
               SCIP_CALL( createExpr(targetset, targetblkmem, &exprcopy, targetexprhdlr, targetexprdata, 0, NULL, ownerdatacreate, ownerdatacreatedata, ownerdatafree) );
            }

            /* let future owner creates its data and store its free callback in the expr */
            SCIP_CALL( createExprOwnerData(targetset, exprcopy, ownerdatacreate, ownerdatacreatedata, ownerdatafree) );

            /* store targetexpr */
            expriteruserdata.ptrval = exprcopy;
            SCIPexpriteratorSetCurrentUserData(it, expriteruserdata);

            break;
         }

         case SCIP_EXPRITER_VISITEDCHILD :
         {
            /* just visited child so a copy of himself should be available; append it */
            SCIP_EXPR* exprcopy;
            SCIP_EXPR* childcopy;

            exprcopy = (SCIP_EXPR*)SCIPexpriteratorGetCurrentUserData(it).ptrval;

            /* get copy of child */
            childcopy = (SCIP_EXPR*)SCIPexpriteratorGetChildUserDataDFS(it).ptrval;
            if( childcopy == NULL )
            {
               /* abort */
               /* release exprcopy (should free also the already copied children) */
               SCIP_CALL( SCIPreleaseExpr(targetscip, (SCIP_EXPR**)&exprcopy) );

               expriteruserdata.ptrval = NULL;
               SCIPexpriteratorSetCurrentUserData(it, expriteruserdata);

               expr = SCIPexpriteratorSkipDFS(it);
               continue;
            }

            /* append child to exprcopy */
            SCIP_CALL( SCIPappendExprChild(targetscip, exprcopy, childcopy) );

            /* release childcopy (still captured by exprcopy) */
            SCIP_CALL( SCIPreleaseExpr(targetscip, &childcopy) );

            break;
         }

         default:
            /* we should never be called in this stage */
            SCIPABORT();
            break;
      }

      expr = SCIPexpriteratorGetNext(it);
   }

   /* the target expression should be stored in the userdata of the sourceexpr (can be NULL if aborted) */
   *targetexpr = (SCIP_EXPR*)SCIPexpriteratorGetExprUserData(it, sourceexpr).ptrval;

   SCIPexpriteratorFree(&it);

   return SCIP_OKAY;
}


/*
 * quadratic representation of expression
 */

/** frees data of quadratic representation of expression, if any */
static
void quadFree(
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_EXPR*            expr                /**< expression whose quadratic data will be released */
   )
{
   int i;
   int n;

   assert(blkmem != NULL);
   assert(expr != NULL);

   expr->quadchecked = FALSE;

   if( expr->quaddata == NULL )
      return;

   n = expr->quaddata->nquadexprs;

   BMSfreeBlockMemoryArrayNull(blkmem, &expr->quaddata->linexprs, expr->quaddata->nlinexprs);
   BMSfreeBlockMemoryArrayNull(blkmem, &expr->quaddata->lincoefs, expr->quaddata->nlinexprs);
   BMSfreeBlockMemoryArrayNull(blkmem, &expr->quaddata->bilinexprterms, expr->quaddata->nbilinexprterms);
   BMSfreeBlockMemoryArrayNull(blkmem, &expr->quaddata->eigenvalues, n);
   BMSfreeBlockMemoryArrayNull(blkmem, &expr->quaddata->eigenvectors, n * n);

   for( i = 0; i < n; ++i )
   {
      BMSfreeBlockMemoryArrayNull(blkmem, &expr->quaddata->quadexprterms[i].adjbilin,
         expr->quaddata->quadexprterms[i].adjbilinsize);
   }
   BMSfreeBlockMemoryArrayNull(blkmem, &expr->quaddata->quadexprterms, n);

   BMSfreeBlockMemory(blkmem, &expr->quaddata);
}

/** first time seen quadratically and
 * seen before linearly --> --nlinterms; assign 2; ++nquadterms
 * not seen before linearly --> assing 1; ++nquadterms
 *
 * seen before --> assign += 1
 */
static
SCIP_RETCODE quadDetectProcessExpr(
   SCIP_EXPR*            expr,               /**< the expression */
   SCIP_HASHMAP*         seenexpr,           /**< hash map */
   int*                  nquadterms,         /**< number of quadratic terms */
   int*                  nlinterms           /**< number of linear terms */
   )
{
   if( SCIPhashmapExists(seenexpr, (void*)expr) )
   {
      int nseen = SCIPhashmapGetImageInt(seenexpr, (void*)expr);

      if( nseen < 0 )
      {
         /* only seen linearly before */
         assert(nseen == -1);

         --*nlinterms;
         ++*nquadterms;
         SCIP_CALL( SCIPhashmapSetImageInt(seenexpr, (void*)expr, 2) );
      }
      else
      {
         assert(nseen > 0);
         SCIP_CALL( SCIPhashmapSetImageInt(seenexpr, (void*)expr, nseen + 1) );
      }
   }
   else
   {
      ++*nquadterms;
      SCIP_CALL( SCIPhashmapInsertInt(seenexpr, (void*)expr, 1) );
   }

   return SCIP_OKAY;
}

/** returns a quadexprterm that contains the expr
 *
 * it either finds one that already exists or creates a new one
 */
static
SCIP_RETCODE quadDetectGetQuadexprterm(
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_EXPR*            expr,               /**< the expression */
   SCIP_HASHMAP*         expr2idx,           /**< map: expr to index in quadexpr->quadexprterms */
   SCIP_HASHMAP*         seenexpr,           /**< map: expr to number of times it was seen */
   SCIP_QUADEXPR*        quadexpr,           /**< data of quadratic representation of expression */
   SCIP_QUADEXPR_QUADTERM** quadexprterm     /**< buffer to store quadexprterm */
   )
{
   assert(expr != NULL);
   assert(expr2idx != NULL);
   assert(quadexpr != NULL);
   assert(quadexprterm != NULL);

   if( SCIPhashmapExists(expr2idx, (void*)expr) )
   {
      *quadexprterm = &quadexpr->quadexprterms[SCIPhashmapGetImageInt(expr2idx, (void*)expr)];
      assert((*quadexprterm)->expr == expr);
   }
   else
   {
      SCIP_CALL( SCIPhashmapInsertInt(expr2idx, expr, quadexpr->nquadexprs) );
      *quadexprterm = &quadexpr->quadexprterms[quadexpr->nquadexprs];
      ++quadexpr->nquadexprs;

      (*quadexprterm)->expr = expr;
      (*quadexprterm)->sqrcoef = 0.0;
      (*quadexprterm)->sqrexpr = NULL;
      (*quadexprterm)->lincoef = 0.0;
      (*quadexprterm)->nadjbilin = 0;
      (*quadexprterm)->adjbilinsize = SCIPhashmapGetImageInt(seenexpr, (void*)expr);
      SCIP_ALLOC( BMSallocBlockMemoryArray(blkmem, &(*quadexprterm)->adjbilin, (*quadexprterm)->adjbilinsize) );
   }

   return SCIP_OKAY;
}

/** @} */


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
 * parse(Expr|Term|Base) returns an SCIP_EXPR
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
SCIP_RETCODE parseExpr(SCIP*, SCIP_CONSHDLR*, SCIP_HASHMAP*, const char*, const char**, SCIP_EXPR**);

/** Parses base to build a value, variable, sum, or function-like ("func(...)") expression.
 * <pre>
 * Base       -> "number" | "<varname>" | "(" Expression ")" | Op "(" OpExpression ")
 * </pre>
 */
static
SCIP_RETCODE parseBase(
   SCIP_SET*             set,                /**< global SCIP settings */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_HASHMAP*         vartoexprvarmap,    /**< hashmap to map between SCIP vars and var expressions */
   const char*           expr,               /**< expr that we are parsing */
   const char**          newpos,             /**< buffer to store the position of expr where we finished reading */
   SCIP_EXPR**           basetree            /**< buffer to store the expr parsed by Base */
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
      SCIP_CALL( SCIPparseVarName(set->scip, expr, &var, (char**)newpos) );

      if( var == NULL )
      {
         SCIPerrorMessage("Could not find variable with name '%s'\n", expr);
         return SCIP_READERROR;
      }
      expr = *newpos;

      /* check if we have already created an expression out of this var */
      if( SCIPhashmapExists(vartoexprvarmap, (void*)var) )
      {
         debugParse("Variable <%s> has already been parsed, capturing its expression\n", SCIPvarGetName(var));
         *basetree = (SCIP_EXPR*)SCIPhashmapGetImage(vartoexprvarmap, (void*)var);
         SCIPcaptureExpr(*basetree);
      }
      else
      {
         debugParse("First time parsing variable <%s>, creating varexpr and adding it to hashmap\n", SCIPvarGetName(var));
         /* intentionally not using createExprVar here, since parsed expressions are not part of a constraint (they will be copied when a constraint is created) */
         SCIP_CALL( SCIPcreateConsExprExprVar(set->scip, basetree, var) );
         SCIP_CALL( SCIPhashmapInsert(vartoexprvarmap, (void*)var, (void*)(*basetree)) );
      }
   }
   else if( *expr == '(' )
   {
      /* parse expression */
      SCIP_CALL( parseExpr(set, blkmem, vartoexprvarmap, ++expr, newpos, basetree) );
      expr = *newpos;

      /* expect ')' */
      if( *expr != ')' )
      {
         SCIPerrorMessage("Read a '(', parsed expression inside --> expecting closing ')'. Got <%c>: rest of string <%s>\n", *expr, expr);
         SCIP_CALL( SCIPreleaseExpr(scip, basetree) );
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
      SCIP_CALL( SCIPcreateConsExprExprValue(scip, basetree, value) );
   }
   else if( isalpha(*expr) )
   {
      /* a (function) name is coming, should find exprhandler with such name */
      int i;
      char operatorname[SCIP_MAXSTRLEN];
      SCIP_EXPRHDLR* exprhdlr;
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

      /* after name we must see a '(' */
      if( *expr != '(' )
      {
         SCIPerrorMessage("Expected '(' after operator name <%s>, but got %s.\n", operatorname, expr);
         return SCIP_READERROR;
      }

      /* search for expression handler */
      exprhdlr = SCIPsetFindExprhdlr(set, operatorname);

      /* check expression handler exists and has a parsing method */
      if( exprhdlr == NULL )
      {
         SCIPerrorMessage("No expression handler with name <%s> found.\n", operatorname);
         return SCIP_READERROR;
      }

      ++expr;
      SCIP_CALL( SCIPcallExprhdlrParse(set->scip, exprhdlr, expr, newpos, basetree, &success) );

      if( !success )
      {
         SCIPerrorMessage("Error while expression handler <%s> was parsing %s\n", operatorname, expr);
         assert(*basetree == NULL);
         return SCIP_READERROR;
      }
      expr = *newpos;

      /* we should see the ')' of Op "(" OpExpression ") */
      assert(*expr == ')');

      /* move one character forward */
      ++expr;
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
   SCIP_SET*             set,                /**< global SCIP settings */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_Bool             isdenominator,      /**< whether factor is in the denominator */
   SCIP_HASHMAP*         vartoexprvarmap,    /**< hashmap to map between scip vars and var expressions */
   const char*           expr,               /**< expr that we are parsing */
   const char**          newpos,             /**< buffer to store the position of expr where we finished reading */
   SCIP_EXPR**           factortree          /**< buffer to store the expr parsed by Factor */
   )
{
   SCIP_EXPR*  basetree;
   SCIP_Real exponent;

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
         SCIP_CALL( SCIPreleaseExpr(scip, &basetree) );
         return SCIP_READERROR;
      }

      if( *expr == '(' )
      {
         ++expr;

         /* it is exponent with parenthesis; expect number possibly starting with + or - */
         if( !SCIPstrToRealValue(expr, &exponent, (char**)&expr) )
         {
            SCIPerrorMessage("error parsing number from <%s>\n", expr);
            SCIP_CALL( SCIPreleaseExpr(scip, &basetree) );
            return SCIP_READERROR;
         }

         /* expect the ')' */
         while( isspace((unsigned char)*expr) )
            ++expr;
         if( *expr != ')' )
         {
            SCIPerrorMessage("error in parsing exponent: expected ')', received <%c> from <%s>\n", *expr,  expr);
            SCIP_CALL( SCIPreleaseExpr(scip, &basetree) );
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
            if( !SCIPstrToRealValue(expr, &exponent, (char**)&expr) )
            {
               SCIPerrorMessage("error parsing number from <%s>\n", expr);
               SCIP_CALL( SCIPreleaseExpr(scip, &basetree) );
               return SCIP_READERROR;
            }
         }
         else
         {
            SCIPerrorMessage("error in parsing exponent, expected a digit, received <%c> from <%s>\n", *expr,  expr);
            SCIP_CALL( SCIPreleaseExpr(scip, &basetree) );
            return SCIP_READERROR;
         }
      }

      debugParse("parsed the exponent %g\n", exponent); /*lint !e506 !e681*/
   }
   else
   {
      /* there is no explicit exponent */
      exponent = 1.0;
   }
   *newpos = expr;

   /* multiply with -1 when we are in the denominator */
   if( isdenominator )
      exponent *= -1.0;

   /* create power */
   if( exponent != 1.0 )
   {
      SCIP_CALL( SCIPcreateConsExprExprPow(set->scip, factortree, basetree, exponent) );
      SCIP_CALL( SCIPreleaseExpr(scip, &basetree) );
   }
   else
      /* Factor consists of this unique Base */
      *factortree = basetree;

   return SCIP_OKAY;
}

/** Parses a term and builds a product-expression, where each factor is a child.
 * <pre>
 * Term -> Factor { ("*" | "/" ) Factor }
 * </pre>
 */
static
SCIP_RETCODE parseTerm(
   SCIP_SET*             set,                /**< global SCIP settings */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_HASHMAP*         vartoexprvarmap,    /**< hashmap to map between scip vars and var expressions */
   const char*           expr,               /**< expr that we are parsing */
   const char**          newpos,             /**< buffer to store the position of expr where we finished reading */
   SCIP_EXPR**           termtree            /**< buffer to store the expr parsed by Term */
   )
{
   SCIP_EXPR* factortree;

   debugParse("parsing term from %s\n", expr);

   /* parse Factor */
   /* ignore whitespace */
   while( isspace((unsigned char)*expr) )
      ++expr;

   SCIP_CALL( parseFactor(set, blkmem, FALSE, vartoexprvarmap, expr, newpos, &factortree) );
   expr = *newpos;

   debugParse("back to parsing Term, continue parsing from %s\n", expr);

   /* check if Terms has another Factor incoming */
   while( isspace((unsigned char)*expr) )
      ++expr;
   if( *expr == '*' || *expr == '/' )
   {
      /* initialize termtree as a product expression with a single term, so we can append the extra Factors */
      SCIP_CALL( SCIPcreateConsExprExprProduct(set->scip, termtree, 1, &factortree, 1.0) );
      SCIP_CALL( SCIPreleaseExpr(scip, &factortree) );

      /* loop: parse Factor, find next symbol */
      do
      {
         SCIP_RETCODE retcode;
         SCIP_Bool isdivision;

         isdivision = (*expr == '/') ? TRUE : FALSE;

         debugParse("while parsing term, read char %c\n", *expr); /*lint !e506 !e681*/

         ++expr;
         retcode = parseFactor(set->scip, isdivision, vartoexprvarmap, expr, newpos, &factortree);

         /* release termtree, if parseFactor fails with a read-error */
         if( retcode == SCIP_READERROR )
         {
            SCIP_CALL( SCIPreleaseExpr(scip, termtree) );
         }
         SCIP_CALL( retcode );

         /* append newly created factor */
         SCIP_CALL( SCIPappendConsExprExprProductExpr(scip, *termtree, factortree) );
         SCIP_CALL( SCIPreleaseExpr(scip, &factortree) );

         /* find next symbol */
         expr = *newpos;
         while( isspace((unsigned char)*expr) )
            ++expr;
      } while( *expr == '*' || *expr == '/' );
   }
   else
   {
      /* Term consists of this unique factor */
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
   SCIP_SET*             set,                /**< global SCIP settings */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_HASHMAP*         vartoexprvarmap,    /**< hashmap to map between scip vars and var expressions */
   const char*           expr,               /**< expr that we are parsing */
   const char**          newpos,             /**< buffer to store the position of expr where we finished reading */
   SCIP_EXPR**           exprtree            /**< buffer to store the expr parsed by Expr */
   )
{
   SCIP_Real sign;
   SCIP_EXPR* termtree;

   debugParse("parsing expression %s\n", expr); /*lint !e506 !e681*/

   /* ignore whitespace */
   while( isspace((unsigned char)*expr) )
      ++expr;

   /* if '+' or '-', store it */
   sign = 1.0;
   if( *expr == '+' || *expr == '-' )
   {
      debugParse("while parsing expression, read char %c\n", *expr); /*lint !e506 !e681*/
      sign = *expr == '+' ? 1.0 : -1.0;
      ++expr;
   }

   SCIP_CALL( parseTerm(set, blkmem, vartoexprvarmap, expr, newpos, &termtree) );
   expr = *newpos;

   debugParse("back to parsing expression (we have the following term), continue parsing from %s\n", expr); /*lint !e506 !e681*/

   /* check if Expr has another Term incoming */
   while( isspace((unsigned char)*expr) )
      ++expr;
   if( *expr == '+' || *expr == '-' )
   {
      if( SCIPisExprValue(termtree) )
      {
         /* initialize exprtree as a sum expression with a constant only, so we can append the following terms */
         SCIP_CALL( SCIPcreateConsExprExprSum(set->scip, exprtree, 0, NULL, NULL, sign * SCIPgetConsExprExprValueValue(termtree)) );
         SCIP_CALL( SCIPreleaseExpr(scip, &termtree) );
      }
      else
      {
         /* initialize exprtree as a sum expression with a single term, so we can append the following terms */
         SCIP_CALL( SCIPcreateConsExprExprSum(set->scip, exprtree, 1, &termtree, &sign, 0.0) );
         SCIP_CALL( SCIPreleaseExpr(scip, &termtree) );
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

         debugParse("while parsing expression, read coefficient %g\n", coef); /*lint !e506 !e681*/

         retcode = parseTerm(scip, conshdlr, vartoexprvarmap, expr, newpos, &termtree);

         /* release exprtree if parseTerm fails with an read-error */
         if( retcode == SCIP_READERROR )
         {
            SCIP_CALL( SCIPreleaseExpr(scip, exprtree) );
         }
         SCIP_CALL( retcode );

         /* append newly created term */
         SCIP_CALL( SCIPappendConsExprExprSumExpr(scip, *exprtree, termtree, coef) );
         SCIP_CALL( SCIPreleaseExpr(scip, &termtree) );

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
         SCIP_CALL( SCIPcreateConsExprExprSum(set->scip, exprtree, 1, &termtree, &sign, 0.0) );
         SCIP_CALL( SCIPreleaseExpr(scip, &termtree) );
      }
      else
         *exprtree = termtree;
   }

   *newpos = expr;

   return SCIP_OKAY;
}

/** @} */  /* end of parsing methods */


/** @name Simplifying expressions (hashing, common subexpressions, simplify)
 * @{

/** returns an equivalent expression for a given expression if possible
 *
 * it adds the expression to key2expr if the map does not contain the key
 */
static
SCIP_RETCODE findEqualExpr(
//   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_EXPR*            expr,               /**< expression to replace */
   SCIP_MULTIHASH*       key2expr,           /**< mapping of hashes to expressions */
   SCIP_EXPR**           newexpr             /**< pointer to store an equivalent expression (NULL if there is none) */
   )
{  /*lint --e{438}*/
   SCIP_MULTIHASHLIST* multihashlist;

//   assert(set != NULL);
   assert(expr != NULL);
   assert(key2expr != NULL);
   assert(newexpr != NULL);

   *newexpr = NULL;
   multihashlist = NULL;
   do
   {
      /* search for an equivalent expression */
      *newexpr = (SCIP_EXPR*)(SCIPmultihashRetrieveNext(key2expr, &multihashlist, (void*)expr));

      if( *newexpr == NULL )
      {
         /* processed all expressions like expr from hash table, so insert expr */
         SCIP_CALL( SCIPmultihashInsert(key2expr, (void*) expr) );
         break;
      }
      else if( expr != *newexpr )
      {
         assert(SCIPexprCompare(expr, *newexpr) == 0);
         break;
      }
      else
      {
         /* can not replace expr since it is already contained in the hashtablelist */
         assert(expr == *newexpr);
         *newexpr = NULL;
         break;
      }
   }
   while( TRUE ); /*lint !e506*/

   return SCIP_OKAY;
}

/** get key of hash element */
static
SCIP_DECL_HASHGETKEY(hashCommonSubexprGetKey)
{
   return elem;
}  /*lint !e715*/

/** checks if two expressions are structurally the same */
static
SCIP_DECL_HASHKEYEQ(hashCommonSubexprEq)
{
   SCIP_EXPR* expr1;
   SCIP_EXPR* expr2;

   expr1 = (SCIP_EXPR*)key1;
   expr2 = (SCIP_EXPR*)key2;
   assert(expr1 != NULL);
   assert(expr2 != NULL);

   return expr1 == expr2 || SCIPexprCompare(expr1, expr2) == 0;
}  /*lint !e715*/

/** get value of hash element when comparing with another expression */
static
SCIP_DECL_HASHKEYVAL(hashCommonSubexprKeyval)
{
   SCIP_EXPR* expr;
   SCIP_EXPRITER* hashiterator;

   expr = (SCIP_EXPR*) key;
   assert(expr != NULL);

   hashiterator = (SCIP_EXPRITER*) userptr;
   assert(hashiterator != NULL);

   return SCIPexpriteratorGetExprUserData(hashiterator, expr).uintval;
}  /*lint !e715*/

/** hashes an expression using an already existing iterator
 *
 * The iterator must by of type DFS with allowrevisit=FALSE and only the leaveexpr stage enabled.
 * The hashes of all visited expressions will be stored in the iterators expression data.
 */
static
SCIP_RETCODE hashExpr(
   SCIP_SET*             set,                /**< global SCIP settings */
   BMS_BUFMEM*           bufmem,             /**< buffer memory */
   SCIP_EXPR*            expr,               /**< expression to hash */
   SCIP_EXPRITER*        hashiterator,       /**< iterator to use for hashing */
   int*                  nvisitedexprs       /**< counter to increment by the number of expressions visited, or NULL */
   )
{
   SCIP_EXPRITER_USERDATA iterdata;
   unsigned int* childrenhashes;
   int childrenhashessize;
   int i;

   assert(set != NULL);
   assert(expr != NULL);
   assert(hashiterator != NULL);

   childrenhashessize = 5;
   SCIP_ALLOC( BMSallocBufferMemoryArray(bufmem, &childrenhashes, childrenhashessize) );

   for( expr = SCIPexpriteratorRestartDFS(hashiterator, expr); !SCIPexpriteratorIsEnd(hashiterator); expr = SCIPexpriteratorGetNext(hashiterator) ) /*lint !e441*/
   {
      assert(SCIPexpriteratorGetStageDFS(hashiterator) == SCIP_EXPRITER_LEAVEEXPR);

      if( nvisitedexprs != NULL )
         ++*nvisitedexprs;

      /* collect hashes of children */
      if( childrenhashessize < expr->nchildren )
      {
         childrenhashessize = SCIPsetCalcMemGrowSize(set, expr->nchildren);
         SCIP_ALLOC( BMSreallocBufferMemoryArray(bufmem, &childrenhashes, childrenhashessize) );
      }
      for( i = 0; i < expr->nchildren; ++i )
         childrenhashes[i] = SCIPexpriteratorGetExprUserData(hashiterator, expr->children[i]).uintval;

      SCIP_CALL( SCIPcallExprhdlrHash(set->scip, expr, &iterdata.uintval, childrenhashes) );

      SCIPexpriteratorSetCurrentUserData(hashiterator, iterdata);
   }

   BMSfreeBufferMemoryArray(bufmem, &childrenhashes);

   return SCIP_OKAY;
}

/** replaces common sub-expressions in a given expression graph by using a hash key for each expression
 *
 *  The algorithm consists of two steps:
 *
 *  1. traverse through all given expressions and compute for each of them a (not necessarily unique) hash
 *
 *  2. initialize an empty hash table and traverse through all expression; check for each of them if we can find a
 *     structural equivalent expression in the hash table; if yes we replace the expression by the expression inside the
 *     hash table, otherwise we add it to the hash table
 *
 *  @note the hash keys of the expressions are used for the hashing inside the hash table; to compute if two expressions
 *  (with the same hash) are structurally the same we use the function SCIPexprCompare()
 */
static
SCIP_RETCODE replaceCommonSubexpressions(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< dynamic problem statistics */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   BMS_BUFMEM*           bufmem,             /**< buffer memory */
   SCIP_EXPR**           exprs,              /**< expressions (possibly replaced by equivalent on output) */
   int                   nexprs,             /**< total number of expressions */
   SCIP_Bool*            replacedroot        /**< buffer to store whether any root expression (expression in exprs) was replaced */
   )
{
   SCIP_EXPRITER* hashiterator;
   SCIP_EXPRITER* repliterator;
   SCIP_MULTIHASH* key2expr;
   SCIP_CONSDATA* consdata;
   int i;
   int nvisitedexprs = 0;

   assert(set != NULL);
   assert(stat != NULL);
   assert(exprs != NULL);
   assert(nexprs >= 0);
   assert(replacedroot != NULL);

   if( nexprs == 0 )
      return SCIP_OKAY;

   SCIP_CALL( SCIPexpriteratorCreate(stat, blkmem, &hashiterator) );
   SCIP_CALL( SCIPexpriteratorInit(hashiterator, NULL, SCIP_EXPRITER_DFS, FALSE) );
   SCIPexpriteratorSetStagesDFS(hashiterator, SCIP_EXPRITER_LEAVEEXPR);

   /* compute all hashes for each sub-expression */
   for( i = 0; i < nexprs; ++i )
   {
      assert(exprs[i] != NULL);
      SCIP_CALL( hashExpr(set, bufmem, exprs[i], hashiterator, &nvisitedexprs) );
   }

   /* replace equivalent sub-expressions */
   SCIP_CALL( SCIPmultihashCreate(&key2expr, blkmem, nvisitedexprs,
         hashCommonSubexprGetKey, hashCommonSubexprEq, hashCommonSubexprKeyval, (void*)hashiterator) );

   SCIP_CALL( SCIPexpriteratorCreate(stat, blkmem, &repliterator) );

   for( i = 0; i < nexprs; ++i )
   {
      SCIP_EXPR* newroot;
      SCIP_EXPR* newchild;
      SCIP_EXPR* child;

      /* check the root for equivalence separately first */
      SCIP_CALL( findEqualExpr(exprs[i], key2expr, &newroot) );

      if( newroot != NULL )
      {
         assert(newroot != exprs[i]);
         assert(SCIPexprCompare(exprs[i], newroot) == 0);

         SCIPdebugMsg(scip, "replacing common root expression of %dth expr: %p -> %p\n", i, (void*)exprs[i], (void*)newroot);

         SCIP_CALL( SCIPreleaseExpr(set->scip, &exprs[i]) );

         exprs[i] = newroot;
         SCIPcaptureExpr(newroot);

         *replacedroot = TRUE;

         continue;
      }

      /* replace equivalent sub-expressions in the tree */
      SCIP_CALL( SCIPexpriteratorInit(repliterator, exprs[i], SCIP_EXPRITER_DFS, FALSE) );
      SCIPexpriteratorSetStagesDFS(repliterator, SCIP_EXPRITER_VISITINGCHILD);

      while( !SCIPexpriteratorIsEnd(repliterator) )
      {
         child = SCIPexpriteratorGetChildExprDFS(repliterator);
         assert(child != NULL);

         /* try to find an equivalent expression */
         SCIP_CALL( findEqualExpr(child, key2expr, &newchild) );

         /* replace child with newchild */
         if( newchild != NULL )
         {
            assert(child != newchild);
            assert(SCIPexprCompare(child, newchild) == 0);

            SCIPdebugMsg(scip, "replacing common child expression %p -> %p\n", (void*)child, (void*)newchild);

            SCIP_CALL( SCIPreplaceExprChild(scip, SCIPexpriteratorGetCurrent(repliterator), SCIPexpriteratorGetChildIdxDFS(repliterator), newchild) );

            (void) SCIPexpriteratorSkipDFS(repliterator);
         }
         else
         {
            (void) SCIPexpriteratorGetNext(repliterator);
         }
      }
   }

   /* free memory */
   SCIPexpriteratorFree(&repliterator);
   SCIPmultihashFree(&key2expr);
   SCIPexpriteratorFree(&hashiterator);

   return SCIP_OKAY;
}

/** helper function to simplify an expression and its subexpressions */
static
SCIP_RETCODE simplifyConsExprExpr(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< dynamic problem statistics */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_EXPR*            rootexpr,           /**< expression to be simplified */
   SCIP_EXPR**           simplified,         /**< buffer to store simplified expression */
   SCIP_Bool*            changed,            /**< buffer to store if rootexpr actually changed */
   SCIP_Bool*            infeasible          /**< buffer to store whether infeasibility has been detected */
   )
{
   SCIP_EXPR* expr;
   SCIP_EXPRITER* it;

   assert(set != NULL);
   assert(stat != NULL);
   assert(blkmem != NULL);
   assert(rootexpr != NULL);
   assert(simplified != NULL);
   assert(changed != NULL);
   assert(infeasible != NULL);

   /* simplify bottom up
    * when leaving an expression it simplifies it and stores the simplified expr in its iterators expression data
    * after the child was visited, it is replaced with the simplified expr
    */
   SCIP_CALL( SCIPexpriteratorCreate(stat, blkmem, &it) );
   SCIP_CALL( SCIPexpriteratorInit(it, rootexpr, SCIP_EXPRITER_DFS, TRUE) );  /* TODO can we set allowrevisited to FALSE?*/
   SCIPexpriteratorSetStagesDFS(it, SCIP_EXPRITER_VISITEDCHILD | SCIP_EXPRITER_LEAVEEXPR);

   *changed = FALSE;
   *infeasible = FALSE;
   for( expr = SCIPexpriteratorGetCurrent(it); !SCIPexpriteratorIsEnd(it); expr = SCIPexpriteratorGetNext(it) ) /*lint !e441*/
   {
      switch( SCIPexpriteratorGetStageDFS(it) )
      {
         case SCIP_EXPRITER_VISITEDCHILD:
         {
            SCIP_EXPR* newchild;
            SCIP_EXPR* child;

            newchild = (SCIP_EXPR*)SCIPexpriteratorGetChildUserDataDFS(it).ptrval;
            child = SCIPexpriteratorGetChildExprDFS(it);
            assert(newchild != NULL);

            /* if child got simplified, replace it with the new child */
            if( newchild != child )
            {
               SCIP_CALL( SCIPreplaceExprChild(scip, expr, SCIPexpriteratorGetChildIdxDFS(it), newchild) );
            }

            /* we do not need to hold newchild anymore */
            SCIP_CALL( SCIPreleaseExpr(scip, &newchild) );

            break;
         }

         case SCIP_EXPRITER_LEAVEEXPR:
         {
            SCIP_EXPR* refexpr = NULL;
            SCIP_EXPRITER_USERDATA iterdata;

            /* TODO we should do constant folding (handle that all children are value-expressions) here in a generic way
             * instead of reimplementing it in every handler
             */

            /* use simplification of expression handlers */
            if( SCIPexprhdlrHasSimplify(expr->exprhdlr) )
            {
               SCIP_CALL( SCIPcallExprhdlrSimplify(scip, conshdlr, expr, &refexpr) );
               if( expr != refexpr )
                  *changed = TRUE;
            }
            else
            {
               /* if an expression handler doesn't implement simplify, we assume all those type of expressions are simplified
                * we have to capture it, since it must simulate a "normal" simplified call in which a new expression is created
                */
               refexpr = expr;
               SCIPcaptureExpr(refexpr);
            }
            assert(refexpr != NULL);

            iterdata.ptrval = (void*) refexpr;
            SCIPexpriteratorSetCurrentUserData(it, iterdata);

            break;
         }

         default:
            SCIPABORT(); /* we should never be called in this stage */
            break;
      }
   }

   *simplified = (SCIP_EXPR*)SCIPexpriteratorGetExprUserData(it, rootexpr).ptrval;
   assert(*simplified != NULL);

   SCIPexpriteratorFree(&it);

   return SCIP_OKAY;
}

/** @} */  /* end of simplify methods */


/** evaluate and forward-differentiate expression */
static
SCIP_RETCODE evalAndDiff(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< dynamic problem statistics */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_EXPR*            expr,               /**< expression to be evaluated */
   SCIP_SOL*             sol,                /**< solution to be evaluated */
   SCIP_Longint          soltag              /**< tag that uniquely identifies the solution (with its values), or 0. */
   )
{
   SCIP_EXPRITER* it;

   assert(set != NULL);
   assert(stat != NULL);
   assert(blkmem != NULL);
   assert(expr != NULL);

   /* assume we'll get a domain error, so we don't have to get this expr back if we abort the iteration
    * if there is no domain error, then we will overwrite the evalvalue in the last leaveexpr stage
    */
   expr->evalvalue = SCIP_INVALID;
   expr->evaltag = soltag;
   expr->dot = SCIP_INVALID;

   SCIP_CALL( SCIPexpriteratorCreate(stat, blkmem, &it) );
   SCIP_CALL( SCIPexpriteratorInit(it, expr, SCIP_EXPRITER_DFS, TRUE) );
   SCIPexpriteratorSetStagesDFS(it, SCIP_EXPRITER_LEAVEEXPR);

   for( expr = SCIPexpriteratorGetCurrent(it); !SCIPexpriteratorIsEnd(it); expr = SCIPexpriteratorGetNext(it) )  /*lint !e441*/
   {
      /* evaluate expression only if necessary */
      if( soltag == 0 || expr->evaltag != soltag )
      {
         SCIP_CALL( SCIPcallExprhdlrEval(set->scip, expr, &expr->evalvalue, NULL, sol) );

         expr->evaltag = soltag;
      }

      if( expr->evalvalue == SCIP_INVALID ) /*lint !e777*/
         break;

      /* compute forward diff */
      SCIP_CALL( SCIPcallExprhdlrFwdiff(set->scip, expr, &expr->dot) );

      if( expr->dot == SCIP_INVALID ) /*lint !e777*/
         break;
   }

   SCIPexpriteratorFree(&it);

   return SCIP_OKAY;
}

/** print statistics for expression handlers */
static
void printExprHdlrStatistics(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler */
   FILE*                 file                /**< file handle, or NULL for standard out */
   )
{
   int i;

   assert(set != NULL);

   SCIPmessageFPrintInfo(messagehdlr, file,
      "Expression Handlers: %10s %10s %10s %10s %10s %10s %10s %10s %10s %10s %10s %10s %10s\n",
      "SimplCalls", "Simplified", "EstimCalls", "#IntEval", "PropCalls", "Cuts", "Cutoffs", "DomReds", "BranchScor", "EstimTime", "PropTime", "IntEvalTi", "SimplifyTi");

   for( i = 0; i < set->nexprhdlrs; ++i )
   {
      SCIP_EXPRHDLR* exprhdlr = set->exprhdlrs[i];
      assert(exprhdlr != NULL);

      SCIPmessageFPrintInfo(messagehdlr, file, "  %-17s:", exprhdlr->name);
      SCIPmessageFPrintInfo(messagehdlr, file, " %10lld", exprhdlr->nsimplifycalls);
      SCIPmessageFPrintInfo(messagehdlr, file, " %10lld", exprhdlr->nsimplified);
      SCIPmessageFPrintInfo(messagehdlr, file, " %10lld", exprhdlr->nestimatecalls);
      SCIPmessageFPrintInfo(messagehdlr, file, " %10lld", exprhdlr->nintevalcalls);
      SCIPmessageFPrintInfo(messagehdlr, file, " %10lld", exprhdlr->npropcalls);
      SCIPmessageFPrintInfo(messagehdlr, file, " %10lld", exprhdlr->ncutsfound);
      SCIPmessageFPrintInfo(messagehdlr, file, " %10lld", exprhdlr->ncutoffs);
      SCIPmessageFPrintInfo(messagehdlr, file, " %10lld", exprhdlr->ndomreds);
      SCIPmessageFPrintInfo(messagehdlr, file, " %10lld", exprhdlr->nbranchscores);
      SCIPmessageFPrintInfo(messagehdlr, file, " %10.2f", SCIPclockGetTime(exprhdlr->estimatetime));
      SCIPmessageFPrintInfo(messagehdlr, file, " %10.2f", SCIPclockGetTime(exprhdlr->proptime));
      SCIPmessageFPrintInfo(messagehdlr, file, " %10.2f", SCIPclockGetTime(exprhdlr->intevaltime));
      SCIPmessageFPrintInfo(messagehdlr, file, " %10.2f", SCIPclockGetTime(exprhdlr->simplifytime));
      SCIPmessageFPrintInfo(messagehdlr, file, "\n");
   }
}



/*
 * Public methods
 */

/**@name Expression Handler Methods */
/**@{ */

/** set the expression handler callbacks to copy and free an expression handler */
void SCIPexprhdlrSetCopyFreeHdlr(
   SCIP_EXPRHDLR*        exprhdlr,           /**< expression handler */
   SCIP_DECL_EXPRCOPYHDLR((*copyhdlr)),      /**< handler copy callback (can be NULL) */
   SCIP_DECL_EXPRFREEHDLR((*freehdlr))       /**< handler free callback (can be NULL) */
)
{
   assert(exprhdlr != NULL);

   exprhdlr->copyhdlr = copyhdlr;
   exprhdlr->freehdlr = freehdlr;
}

/** set the expression handler callbacks to copy and free expression data */
void SCIPexprhdlrSetCopyFreeData(
   SCIP_EXPRHDLR*        exprhdlr,           /**< expression handler */
   SCIP_DECL_EXPRCOPYDATA((*copydata)),      /**< expression data copy callback (can be NULL for expressions without data) */
   SCIP_DECL_EXPRFREEDATA((*freedata))       /**< expression data free callback (can be NULL if data does not need to be freed) */
)
{  /*lint --e{715}*/
   assert(exprhdlr != NULL);

   exprhdlr->copydata = copydata;
   exprhdlr->freedata = freedata;
}

/** set the print callback of an expression handler */
void SCIPexprhdlrSetPrint(
   SCIP_EXPRHDLR*        exprhdlr,           /**< expression handler */
   SCIP_DECL_EXPRPRINT((*print))             /**< print callback (can be NULL) */
)
{
   assert(exprhdlr != NULL);

   exprhdlr->print = print;
}

/** set the parse callback of an expression handler */
void SCIPexprhdlrSetParse(
   SCIP_EXPRHDLR*        exprhdlr,           /**< expression handler */
   SCIP_DECL_EXPRPARSE((*parse))             /**< parse callback (can be NULL) */
)
{
   assert(exprhdlr != NULL);

   exprhdlr->parse = parse;
}

/** set the curvature detection callback of an expression handler */
void SCIPexprhdlrSetCurvature(
   SCIP_EXPRHDLR*        exprhdlr,           /**< expression handler */
   SCIP_DECL_EXPRCURVATURE((*curvature))     /**< curvature detection callback (can be NULL) */
)
{
   assert(exprhdlr != NULL);

   exprhdlr->curvature = curvature;
}

/** set the monotonicity detection callback of an expression handler */
void SCIPexprhdlrSetMonotonicity(
   SCIP_EXPRHDLR*        exprhdlr,           /**< expression handler */
   SCIP_DECL_EXPRMONOTONICITY((*monotonicity)) /**< monotonicity detection callback (can be NULL) */
)
{
   assert(exprhdlr != NULL);

   exprhdlr->monotonicity = monotonicity;
}

/** set the integrality detection callback of an expression handler */
void SCIPexprhdlrSetIntegrality(
   SCIP_EXPRHDLR*        exprhdlr,           /**< expression handler */
   SCIP_DECL_EXPRINTEGRALITY((*integrality)) /**< integrality detection callback (can be NULL) */
)
{
   assert(exprhdlr != NULL);

   exprhdlr->integrality = integrality;
}

/** set the hash callback of an expression handler */
void SCIPexprhdlrSetHash(
   SCIP_EXPRHDLR*        exprhdlr,           /**< expression handler */
   SCIP_DECL_EXPRHASH((*hash))               /**< hash callback (can be NULL) */
)
{
   assert(exprhdlr != NULL);

   exprhdlr->hash = hash;
}

/** set the compare callback of an expression handler */
void SCIPexprhdlrSetCompare(
   SCIP_EXPRHDLR*        exprhdlr,           /**< expression handler */
   SCIP_DECL_EXPRCOMPARE((*compare))         /**< compare callback (can be NULL) */
)
{
   assert(exprhdlr != NULL);

   exprhdlr->compare = compare;
}

/** set derivative evaluation callbacks of an expression handler */
void SCIPexprhdlrSetDiff(
   SCIP_EXPRHDLR*        exprhdlr,           /**< expression handler */
   SCIP_DECL_EXPRBWDIFF((*bwdiff)),          /**< backward derivative evaluation callback (can be NULL) */
   SCIP_DECL_EXPRFWDIFF((*fwdiff)),          /**< forward derivative evaluation callback (can be NULL) */
   SCIP_DECL_EXPRBWFWDIFF((*bwfwdiff))       /**< backward-forward derivative evaluation callback (can be NULL) */
)
{
   assert(exprhdlr != NULL);

   exprhdlr->bwdiff = bwdiff;
   exprhdlr->fwdiff = fwdiff;
   exprhdlr->bwfwdiff = bwfwdiff;
}

/** set the interval evaluation callback of an expression handler */
void SCIPexprhdlrSetIntEval(
   SCIP_EXPRHDLR*        exprhdlr,           /**< expression handler */
   SCIP_DECL_EXPRINTEVAL((*inteval))         /**< interval evaluation callback (can be NULL) */
)
{
   assert(exprhdlr != NULL);

   exprhdlr->inteval = inteval;
}

/** set the simplify callback of an expression handler */
void SCIPexprhdlrSetSimplify(
   SCIP_EXPRHDLR*        exprhdlr,           /**< expression handler */
   SCIP_DECL_EXPRSIMPLIFY((*simplify))       /**< simplify callback (can be NULL) */
)
{
   assert(exprhdlr != NULL);

   exprhdlr->simplify = simplify;
}

/** set the reverse propagation callback of an expression handler */
void SCIPexprhdlrSetReverseProp(
   SCIP_EXPRHDLR*        exprhdlr,           /**< expression handler */
   SCIP_DECL_EXPRREVERSEPROP((*reverseprop)) /**< reverse propagation callback (can be NULL) */
)
{
   assert(exprhdlr != NULL);

   exprhdlr->reverseprop = reverseprop;
}

/** set the and estimation callbacks of an expression handler */
void SCIPexprhdlrSetEstimate(
   SCIP_EXPRHDLR*        exprhdlr,                /**< expression handler */
   SCIP_DECL_EXPRINITESTIMATES((*initestimates)), /**< initial estimators callback (can be NULL) */
   SCIP_DECL_EXPRESTIMATE((*estimate))            /**< estimator callback (can be NULL) */
)
{
   assert(exprhdlr != NULL);

   exprhdlr->initestimates = initestimates;
   exprhdlr->estimate = estimate;
}

/** gives the name of an expression handler */
const char* SCIPexprhdlrGetName(
   SCIP_EXPRHDLR*        exprhdlr            /**< expression handler */
)
{
   assert(exprhdlr != NULL);

   return exprhdlr->name;
}

/** gives the description of an expression handler (can be NULL) */
const char* SCIPexprhdlrGetDescription(
   SCIP_EXPRHDLR*        exprhdlr            /**< expression handler */
)
{
   assert(exprhdlr != NULL);

   return exprhdlr->desc;
}

/** gives the precedence of an expression handler */
unsigned int SCIPexprhdlrGetPrecedence(
   SCIP_EXPRHDLR*        exprhdlr            /**< expression handler */
)
{
   assert(exprhdlr != NULL);

   return exprhdlr->precedence;
}

/** gives the data of an expression handler */
SCIP_EXPRHDLRDATA* SCIPexprhdlrGetData(
   SCIP_EXPRHDLR*        exprhdlr            /**< expression handler */
)
{
   assert(exprhdlr != NULL);

   return exprhdlr->data;
}

/** returns whether expression handler implements the print callback */
SCIP_Bool SCIPexprhdlrHasPrint(
   SCIP_EXPRHDLR*        exprhdlr            /**< expression handler */
)
{
   assert(exprhdlr != NULL);

   return exprhdlr->print != NULL;
}

/** returns whether expression handler implements the backward differentiation callback */
SCIP_Bool SCIPexprhdlrHasBwdiff(
   SCIP_EXPRHDLR*        exprhdlr            /**< expression handler */
)
{
   assert(exprhdlr != NULL);

   return exprhdlr->bwdiff != NULL;
}

/** returns whether expression handler implements the interval evaluation callback */
SCIP_Bool SCIPexprhdlrHasIntEval(
   SCIP_EXPRHDLR*        exprhdlr            /**< expression handler */
   )
{
   assert(exprhdlr != NULL);

   return exprhdlr->inteval != NULL;
}

/** returns whether expression handler implements the estimator callback */
SCIP_Bool SCIPexprhdlrHasEstimate(
   SCIP_EXPRHDLR*        exprhdlr            /**< expression handler */
   )
{
   assert(exprhdlr != NULL);

   return exprhdlr->estimate != NULL;
}

/** returns whether expression handler implements the initial estimators callback */
SCIP_Bool SCIPexprhdlrHasInitEstimates(
   SCIP_EXPRHDLR*        exprhdlr            /**< expression handler */
   )
{
   assert(exprhdlr != NULL);

   return exprhdlr->initestimates != NULL;
}

/** returns whether expression handler implements the simplification callback */
SCIP_Bool SCIPexprhdlrHasSimplify(
   SCIP_EXPRHDLR*        exprhdlr            /**< expression handler */
   )
{
   assert(exprhdlr != NULL);

   return exprhdlr->simplify != NULL;
}

/** returns whether expression handler implements the curvature callback */
SCIP_Bool SCIPexprhdlrHasCurvature(
   SCIP_EXPRHDLR*        exprhdlr            /**< expression handler */
   )
{
   assert(exprhdlr != NULL);

   return exprhdlr->curvature != NULL;
}

/** returns whether expression handler implements the reverse propagation callback */
SCIP_Bool SCIPexprhdlrHasReverseProp(
   SCIP_EXPRHDLR*        exprhdlr            /**< expression handler */
   )
{
   assert(exprhdlr != NULL);

   return exprhdlr->reverseprop != NULL;
}

/** increments the branching score count of an expression handler */
void SCIPexprhdlrIncrementNBranchScore(
   SCIP_EXPRHDLR*        exprhdlr            /**< expression handler */
   )
{
   assert(exprhdlr != NULL);

   ++exprhdlr->nbranchscores;
}


/** copies the given expression handler to a new scip */
SCIP_RETCODE SCIPexprhdlrCopyInclude(
   SCIP_EXPRHDLR*        exprhdlr,           /**< expression handler */
   SCIP_SET*             targetset           /**< SCIP_SET of SCIP to copy to */
   )
{
   assert(exprhdlr != NULL);
   assert(targetset != NULL);
   assert(targetset->scip != NULL);

   if( exprhdlr->copyhdlr != NULL )
   {
      SCIPsetDebugMsg(set, "including expression handler <%s> in subscip %p\n", SCIPexprhdlrGetName(exprhdlr), (void*)targetset->scip);
      SCIP_CALL( exprhdlr->copyhdlr(targetset->scip, exprhdlr) );
   }
   else
   {
      SCIPsetDebugMsg(set, "expression handler <%s> cannot be copied to subscip %p due to missing copyhdlr callback\n", SCIPexprhdlrGetName(exprhdlr), (void*)targetset->scip);
   }

   return SCIP_OKAY;
}

/** calls the print callback of an expression handler
 *
 * the method prints an expression
 * it is called while iterating over the expression graph at different stages
 */
SCIP_RETCODE SCIPexprhdlrPrintExpr(
   SCIP_EXPRHDLR*        exprhdlr,           /**< expression handler */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler */
   SCIP_EXPR*            expr,               /**< expression */
   SCIP_EXPRITER_STAGE   stage,              /**< stage of expression iteration */
   int                   currentchild,       /**< index of current child if in stage visitingchild or visitedchild */
   unsigned int          parentprecedence,   /**< precedence of parent */
   FILE*                 file                /**< the file to print to */
   )
{
   assert(exprhdlr != NULL);
   assert(set != NULL);
   assert(expr != NULL);
   assert(expr->exprhdlr == exprhdlr);
   assert(messagehdlr != NULL);
   assert(currentchild >= 0);
   assert(currentchild < expr->nchildren);

   if( SCIPexprhdlrHasPrint(exprhdlr) )
   {
      SCIP_CALL( exprhdlr->print(set->scip, expr, stage, currentchild, parentprecedence, file) );
   }
   else
   {
      /* default: <hdlrname>(<child1>, <child2>, ...) */
      switch( stage )
      {
         case SCIP_EXPRITER_ENTEREXPR :
         {
            SCIPmessageFPrintInfo(messagehdlr, file, SCIPexprhdlrGetName(expr->exprhdlr));
            if( expr->nchildren > 0 )
            {
               SCIPmessageFPrintInfo(messagehdlr, file, "(");
            }
            break;
         }

         case SCIP_EXPRITER_VISITEDCHILD :
         {
            if( currentchild < expr->nchildren-1 )
            {
               SCIPmessageFPrintInfo(messagehdlr, file, ", ");
            }
            else
            {
               SCIPmessageFPrintInfo(messagehdlr, file, ")");
            }

            break;
         }

         case SCIP_EXPRITER_VISITINGCHILD :
         case SCIP_EXPRITER_LEAVEEXPR :
         default:
            break;
      }
   }

   return SCIP_OKAY;
}

/** calls the parse callback of an expression handler
 *
 * The method parses an expression.
 * It should be called when parsing an expression and an operator with the expr handler name is found.
 */
SCIP_RETCODE SCIPexprhdlrParseExpr(
   SCIP_EXPRHDLR*        exprhdlr,           /**< expression handler */
   SCIP_SET*             set,                /**< global SCIP settings */
   const char*           string,             /**< string containing expression to be parse */
   const char**          endstring,          /**< buffer to store the position of string after parsing */
   SCIP_EXPR**           expr,               /**< buffer to store the parsed expression */
   SCIP_Bool*            success             /**< buffer to store whether the parsing was successful or not */
   )
{
   assert(exprhdlr != NULL);
   assert(set != NULL);
   assert(expr != NULL);

   *expr = NULL;

   if( exprhdlr->parse == NULL )
   {
      /* TODO we could just look for a comma separated list of operands and try to initialize the expr with this one?
       * That would be sufficient for sin, cos, exp, log, abs, for example.
       */
      SCIPdebugMessage("Expression handler <%s> has no parsing method.\n", SCIPexprhdlrGetName(exprhdlr));
      *success = FALSE;
      return SCIP_OKAY;
   }

   /* give control to exprhdlr's parser */
   SCIP_CALL( exprhdlr->parse(set->scip, exprhdlr, string, endstring, expr, success) );

   assert(*success || (*expr == NULL));

   return SCIP_OKAY;
}

/** calls the curvature check callback of an expression handler
 *
 * See @ref SCIP_DECL_EXPRCURVATURE for details.
 */
SCIP_RETCODE SCIPexprhdlrCurvatureExpr(
   SCIP_EXPRHDLR*        exprhdlr,           /**< expression handler */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_EXPR*            expr,               /**< expression to check the curvature for */
   SCIP_EXPRCURV         exprcurvature,      /**< desired curvature of this expression */
   SCIP_Bool*            success,            /**< buffer to store whether the desired curvature be obtained */
   SCIP_EXPRCURV*        childcurv           /**< array to store required curvature for each child */
   )
{
   assert(exprhdlr != NULL);
   assert(set != NULL);
   assert(expr != NULL);
   assert(expr->exprhdlr == exprhdlr);
   assert(success != NULL);

   *success = FALSE;

   if( exprhdlr->curvature != NULL )
   {
      SCIP_CALL( exprhdlr->curvature(set->scip, expr, exprcurvature, success, childcurv) );
   }

   return SCIP_OKAY;
}

/** calls the hash callback of an expression handler
 *
 * The method hashes an expression by taking the hashes of its children into account.
 */
SCIP_RETCODE SCIPexprhdlrHashExpr(
   SCIP_EXPRHDLR*        exprhdlr,           /**< expression handler */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_EXPR*            expr,               /**< expression to be hashed */
   unsigned int*         hashkey,            /**< buffer to store the hash value */
   unsigned int*         childrenhashes      /**< array with hash values of children */
   )
{
   assert(exprhdlr != NULL);
   assert(set != NULL);
   assert(expr != NULL);
   assert(expr->exprhdlr == exprhdlr);
   assert(hashkey != NULL);
   assert(childrenhashes != NULL || expr->nchildren == 0);

   if( expr->exprhdlr->hash != NULL )
   {
      SCIP_CALL( expr->exprhdlr->hash(set->scip, expr, hashkey, childrenhashes) );
   }
   else
   {
      int i;

      /* compute initial hash from expression handler name if callback is not implemented
       * this can lead to more collisions and thus a larger number of expensive expression compare calls
       */
      *hashkey = 0;
      for( i = 0; expr->exprhdlr->name[i] != '\0'; i++ )
         *hashkey += (unsigned int) expr->exprhdlr->name[i]; /*lint !e571*/

      *hashkey = SCIPcalcFibHash((SCIP_Real)*hashkey);

      /* now make use of the hashkeys of the children */
      for( i = 0; i < expr->nchildren; ++i )
         *hashkey ^= childrenhashes[i];
   }

   return SCIP_OKAY;
}

/** calls the compare callback of an expression handler
 *
 * The method receives two expressions, expr1 and expr2, and returns
 * - -1 if expr1 < expr2
 * - 0  if expr1 = expr2
 * - 1  if expr1 > expr2
 */
int SCIPexprhdlrCompareExpr(
   SCIP_EXPR*            expr1,              /**< first expression in comparison */
   SCIP_EXPR*            expr2               /**< second expression in comparison */
   )
{
   int i;

   assert(expr1 != NULL);
   assert(expr2 != NULL);
   assert(expr1->exprhdlr == expr2->exprhdlr);

   if( expr1->exprhdlr->compare != NULL )
   {
      /* enforces OR1-OR4 */
      return expr1->exprhdlr->compare(expr1, expr2);
   }

   /* enforces OR5: default comparison method of expressions of the same type:
    * expr1 < expr2 if and only if expr1_i = expr2_i for all i < k and expr1_k < expr2_k.
    * if there is no such k, use number of children to decide
    * if number of children is equal, both expressions are equal
    * @note: Warning, this method doesn't know about expression data. So if your expressions have special data,
    * you must implement the compare callback: SCIP_DECL_EXPRCOMPARE
    */
   for( i = 0; i < expr1->nchildren && i < expr2->nchildren; ++i )
   {
      int compareresult = SCIPexprCompare(expr1->children[i], expr2->children[i]);
      if( compareresult != 0 )
         return compareresult;
   }

   return expr1->nchildren == expr2->nchildren ? 0 : expr1->nchildren < expr2->nchildren ? -1 : 1;
}

/** calls the evaluation callback of an expression handler
 *
 * The method evaluates an expression by taking the values of its children into account.
 *
 * Further, allows to evaluate w.r.t. given expression and children values instead of those stored in children expressions.
 */
SCIP_RETCODE SCIPexprhdlrEvalExpr(
   SCIP_EXPRHDLR*        exprhdlr,           /**< expression handler */
   SCIP_SET*             set,                /**< global SCIP settings */
   BMS_BUFMEM*           bufmem,             /**< buffer memory, can be NULL if childrenvals is NULL */
   SCIP_EXPR*            expr,               /**< expression to be evaluated */
   SCIP_Real*            val,                /**< buffer to store value of expression */
   SCIP_Real*            childrenvals,       /**< values for children, or NULL if values stored in children should be used */
   SCIP_SOL*             sol                 /**< solution that is evaluated (can be NULL) */
)
{
   SCIP_Real* origvals = NULL;

   assert(exprhdlr != NULL);
   assert(set != NULL);
   assert(expr != NULL);
   assert(expr->exprhdlr == exprhdlr);
   assert(exprhdlr->eval != NULL);
   assert(val != NULL);

   /* temporarily overwrite the evalvalue in all children with values from childrenvals */
   if( childrenvals != NULL && expr->nchildren > 0 )
   {
      int c;

      assert(bufmem != NULL);

      SCIP_ALLOC( BMSallocBufferMemoryArray(bufmem, &origvals, expr->nchildren) );

      for( c = 0; c < expr->nchildren; ++c )
      {
         origvals[c] = expr->children[c]->evalvalue;
         expr->children[c]->evalvalue = childrenvals[c];
      }
   }

   /* call expression eval callback */
   SCIP_CALL( exprhdlr->eval(set->scip, expr, val, sol) );

   /* if there was some evaluation error (e.g., overflow) that hasn't been caught yet, then do so now */
   if( !SCIPisFinite(*val) )
      *val = SCIP_INVALID;

   /* restore original evalvalues in children */
   if( origvals != NULL )
   {
      int c;
      for( c = 0; c < expr->nchildren; ++c )
         expr->children[c]->evalvalue = origvals[c];

      BMSfreeBufferMemoryArray(bufmem, &origvals);
   }

   return SCIP_OKAY;
}

/** calls the backward derivative evaluation callback of an expression handler
 *
 * The method should compute the partial derivative of expr w.r.t its child at childidx.
 * That is, it returns
 * \f[
 *   \frac{\partial \text{expr}}{\partial \text{child}_{\text{childidx}}}
 * \f]
 *
 * Further, allows to differentiate w.r.t. given expression and children values instead of those stored in children expressions.
 */
SCIP_RETCODE SCIPexprhdlrBwDiffExpr(
   SCIP_EXPRHDLR*        exprhdlr,           /**< expression handler */
   SCIP_SET*             set,                /**< global SCIP settings */
   BMS_BUFMEM*           bufmem,             /**< buffer memory, can be NULL if childrenvals is NULL */
   SCIP_EXPR*            expr,               /**< expression to be differentiated */
   int                   childidx,           /**< index of the child */
   SCIP_Real*            derivative,         /**< buffer to store the partial derivative w.r.t. the i-th children */
   SCIP_Real*            childrenvals,       /**< values for children, or NULL if values stored in children should be used */
   SCIP_Real             exprval             /**< value for expression, used only if childrenvals is not NULL */
   )
{
   SCIP_Real* origchildrenvals;
   SCIP_Real origexprval;
   int c;

   assert(exprhdlr != NULL);
   assert(set != NULL);
   assert(expr != NULL);
   assert(expr->exprhdlr == exprhdlr);
   assert(derivative != NULL);

   if( exprhdlr->bwdiff == NULL )
   {
      *derivative = SCIP_INVALID;
      return SCIP_OKAY;
   }

   if( childrenvals != NULL )
   {
      /* temporarily overwrite the evalvalue in all children and expr with values from childrenvals and exprval, resp. */
      if( expr->nchildren > 0 )
      {
         assert(bufmem != NULL);
         SCIP_ALLOC( BMSallocBufferMemoryArray(bufmem, &origchildrenvals, expr->nchildren) );

         for( c = 0; c < expr->nchildren; ++c )
         {
            origchildrenvals[c] = expr->children[c]->evalvalue;
            expr->children[c]->evalvalue = childrenvals[c];
         }
      }

      origexprval = expr->evalvalue;
      expr->evalvalue = exprval;
   }

   SCIP_CALL( expr->exprhdlr->bwdiff(set->scip, expr, childidx, derivative) );

   /* if there was some evaluation error (e.g., overflow) that hasn't been caught yet, then do so now */
   if( !SCIPisFinite(*derivative) )
      *derivative = SCIP_INVALID;

   /* restore original evalvalues in children */
   if( childrenvals != NULL )
   {
      if( expr->nchildren > 0 )
      {
         for( c = 0; c < expr->nchildren; ++c )
            expr->children[c]->evalvalue = origchildrenvals[c];  /*lint !e644*/

         BMSfreeBufferMemoryArray(bufmem, &origchildrenvals);
      }

      expr->evalvalue = origexprval;   /*lint !e644*/
   }

   return SCIP_OKAY;
}

/** calls the forward differentiation callback of an expression handler
 *
 * See @ref SCIP_DECL_EXPRFWDIFF for details.
 */
SCIP_RETCODE SCIPexprhdlrFwDiffExpr(
   SCIP_EXPRHDLR*        exprhdlr,           /**< expression handler */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_EXPR*            expr,               /**< expression to be differentiated */
   SCIP_Real*            dot,                /**< buffer to store derivative value */
   SCIP_SOL*             direction           /**< direction of the derivative (useful only for var expressions) */
   )
{
   assert(exprhdlr != NULL);
   assert(set != NULL);
   assert(expr != NULL);
   assert(expr->exprhdlr == exprhdlr);
   assert(dot != NULL);

   if( exprhdlr->fwdiff == NULL )
   {
      *dot = SCIP_INVALID;
      return SCIP_OKAY;
   }

   SCIP_CALL( exprhdlr->fwdiff(set->scip, expr, dot, direction) );

   /* if there was some evaluation error (e.g., overflow) that hasn't been caught yet, then do so now */
   if( !SCIPisFinite(*dot) )
      *dot = SCIP_INVALID;

   return SCIP_OKAY;
}

/** calls the evaluation callback for Hessian directions (backward over forward) of an expression handler
 *
 * See @ref SCIP_DECL_EXPRBWFWDIFF for details.
 */
SCIP_RETCODE SCIPexprhdlrBwFwDiffExpr(
   SCIP_EXPRHDLR*        exprhdlr,           /**< expression handler */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_EXPR*            expr,               /**< expression to be differentiated */
   int                   childidx,           /**< index of the child */
   SCIP_Real*            bardot,             /**< buffer to store derivative value */
   SCIP_SOL*             direction           /**< direction of the derivative (useful only for var expressions) */
   )
{
   assert(exprhdlr != NULL);
   assert(set != NULL);
   assert(expr != NULL);
   assert(expr->exprhdlr == exprhdlr);
   assert(childidx >= 0);
   assert(childidx < expr->nchildren);
   assert(bardot != NULL);

   if( exprhdlr->bwfwdiff == NULL )
   {
      *bardot = SCIP_INVALID;
      return SCIP_OKAY;
   }

   SCIP_CALL( expr->exprhdlr->bwfwdiff(set->scip, expr, childidx, bardot, direction) );

   /* if there was some evaluation error (e.g., overflow) that hasn't been caught yet, then do so now */
   if( !SCIPisFinite(*bardot) )
      *bardot = SCIP_INVALID;

   return SCIP_OKAY;
}

/** calls the interval evaluation callback of an expression handler */
SCIP_RETCODE SCIPexprhdlrIntEvalExpr(
   SCIP_EXPRHDLR*        exprhdlr,           /**< expression handler */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_EXPR*            expr,               /**< expression to be evaluated */
   SCIP_INTERVAL*        interval,           /**< buffer where to store interval */
   SCIP_DECL_EXPR_INTEVALVAR((*intevalvar)), /**< callback to be called when interval-evaluating a variable */
   void*                 intevalvardata      /**< data to be passed to intevalvar callback */
   )
{
   assert(exprhdlr != NULL);
   assert(set != NULL);
   assert(expr != NULL);
   assert(expr->exprhdlr == exprhdlr);
   assert(interval != NULL);

   if( exprhdlr->inteval != NULL )
   {
      SCIP_CALL( SCIPclockStart(exprhdlr->intevaltime, set) );
      SCIP_CALL( exprhdlr->inteval(set->scip, expr, interval, intevalvar, intevalvardata) );
      SCIP_CALL( SCIPclockStop(exprhdlr->intevaltime, set) );

      ++exprhdlr->nintevalcalls;
   }

   return SCIP_OKAY;
}

/** calls the estimator callback of an expression handler
 *
 * See @ref SCIP_DECL_EXPRESTIMATE for details.
 */
SCIP_RETCODE SCIPexprhdlrEstimateExpr(
   SCIP_EXPRHDLR*        exprhdlr,           /**< expression handler */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_EXPR*            expr,               /**< expression to be estimated */
   SCIP_SOL*             sol,                /**< solution at which to estimate (NULL for the LP solution) */
   SCIP_Bool             overestimate,       /**< whether the expression needs to be over- or underestimated */
   SCIP_Real             targetvalue,        /**< a value that the estimator shall exceed, can be +/-infinity */
   SCIP_Real*            coefs,              /**< array to store coefficients of estimator */
   SCIP_Real*            constant,           /**< buffer to store constant part of estimator */
   SCIP_Bool*            islocal,            /**< buffer to store whether estimator is valid locally only */
   SCIP_Bool*            success,            /**< buffer to indicate whether an estimator could be computed */
   SCIP_Bool*            branchcand          /**< array to indicate which children (not) to consider for branching */
   )
{
   assert(exprhdlr != NULL);
   assert(set != NULL);
   assert(expr != NULL);
   assert(expr->exprhdlr == exprhdlr);
   assert(coefs != NULL);
   assert(islocal != NULL);
   assert(success != NULL);

   *success = FALSE;

   if( exprhdlr->estimate != NULL )
   {
      SCIP_CALL( SCIPclockStart(exprhdlr->estimatetime, set) );
      SCIP_CALL( exprhdlr->estimate(set->scip, expr, sol, overestimate, targetvalue, coefs, constant, islocal, success, branchcand) );
      SCIP_CALL( SCIPclockStop(exprhdlr->estimatetime, set) );

      /* update statistics */
      ++exprhdlr->nestimatecalls;
   }

   return SCIP_OKAY;
}

/** calls the intitial estimators callback of an expression handler
 *
 * See @ref SCIP_DECL_EXPRINITESTIMATES for details.
 */
SCIP_RETCODE SCIPexprhdlrInitEstimatesExpr(
   SCIP_EXPRHDLR*        exprhdlr,           /**< expression handler */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_EXPR*            expr,               /**< expression to be estimated */
   SCIP_Bool             overestimate,       /**< whether the expression shall be overestimated or underestimated */
   SCIP_Real*            coefs[SCIP_EXPR_MAXINITESTIMATES], /**< buffer to store coefficients of computed estimators */
   SCIP_Real*            constant[SCIP_EXPR_MAXINITESTIMATES], /**< buffer to store constant of computed estimators */
   SCIP_Bool*            islocal[SCIP_EXPR_MAXINITESTIMATES], /**< buffer to return whether estimator validity depends on children activity */
   int*                  nreturned           /**< buffer to store number of estimators that have been computed */
   )
{
   assert(exprhdlr != NULL);
   assert(set != NULL);
   assert(expr != NULL);
   assert(expr->exprhdlr == exprhdlr);
   assert(nreturned != NULL);

   *nreturned = 0;

   if( exprhdlr->initestimates )
   {
      SCIP_CALL( SCIPclockStart(expr->exprhdlr->estimatetime, set) );
      SCIP_CALL( exprhdlr->initestimates(set->scip, expr, overestimate, coefs, constant, islocal, nreturned) );
      SCIP_CALL( SCIPclockStop(expr->exprhdlr->estimatetime, set) );

      ++exprhdlr->nestimatecalls;
   }

   return SCIP_OKAY;
}

/** calls the simplification callback of an expression handler
 *
 * The function receives the expression to be simplified and a pointer to store the simplified expression.
 */
SCIP_RETCODE SCIPexprhdlrSimplifyExpr(
   SCIP_EXPRHDLR*        exprhdlr,           /**< expression handler */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_EXPR*            expr,               /**< expression to simplify */
   SCIP_EXPR**           simplifiedexpr      /**< buffer to store the simplified expression */
   )
{
   assert(exprhdlr != NULL);
   assert(set != NULL);
   assert(expr != NULL);
   assert(expr->exprhdlr == exprhdlr);
   assert(simplifiedexpr != NULL);

   if( exprhdlr->simplify != NULL )
   {
      SCIP_CALL( SCIPclockStart(expr->exprhdlr->simplifytime, set) );
      SCIP_CALL( exprhdlr->simplify(set->scip, expr, simplifiedexpr) );
      SCIP_CALL( SCIPclockStop(expr->exprhdlr->simplifytime) );

      /* update statistics */
      ++exprhdlr->nsimplifycalls;
      if( expr != *simplifiedexpr )
         ++exprhdlr->nsimplified;
   }

   return SCIP_OKAY;
}

/** calls the reverse propagation callback of an expression handler
 *
 * The method propagates given bounds over the children of an expression.
 */
SCIP_RETCODE SCIPexprhdlrReversePropExpr(
   SCIP_EXPRHDLR*        exprhdlr,           /**< expression handler */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_EXPR*            expr,               /**< expression to propagate */
   SCIP_INTERVAL         bounds,             /**< the bounds on the expression that should be propagated */
   SCIP_INTERVAL*        childrenbounds,     /**< array to store computed bounds for children, initialized with current activity */
   SCIP_Bool*            infeasible          /**< buffer to store whether a children bounds were propagated to an empty interval */
   )
{
   assert(exprhdlr != NULL);
   assert(set != NULL);
   assert(expr != NULL);
   assert(expr->exprhdlr == exprhdlr);
   assert(childrenbounds != NULL || expr->nchildren == 0);
   assert(infeasible != NULL);

   *infeasible = FALSE;

   if( exprhdlr->reverseprop != NULL )
   {
      SCIP_CALL( SCIPclockStart(exprhdlr->proptime, set) );
      SCIP_CALL( exprhdlr->reverseprop(set->scip, expr, bounds, childrenbounds, infeasible) );
      SCIP_CALL( SCIPclockStop(exprhdlr->proptime, set) );

      /* update statistics */
      if( *infeasible )
         ++expr->exprhdlr->ncutoffs;
      ++expr->exprhdlr->npropcalls;

      //TODO move into nlhldr_default:
      //assert(*nreductions >= 0);
      //exprhdlr->ndomreds += *nreductions;
   }

   return SCIP_OKAY;
}

/**@} */
