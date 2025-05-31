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

/**@file   objexprhdlr.cpp
 * @brief  C++ wrapper for expression handlers
 * @author Kevin Kofler
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <cassert>

#include "objexprhdlr.h"




/*
 * Data structures
 */

/** expression handler data */
struct SCIP_ExprhdlrData
{
   scip::ObjExprhdlr*    objexprhdlr;        /**< expression handler object */
   SCIP_Bool             deleteobject;       /**< should the expression handler object be deleted when exprhdlr is freed? */
};




/*
 * Callback methods of expression handler
 */

extern "C"
{

/** copy method for expression handler plugins (called when SCIP copies plugins) */
static
SCIP_DECL_EXPRCOPYHDLR(exprCopyhdlrObj)
{  /*lint --e{715}*/
   SCIP_EXPRHDLRDATA* exprhdlrdata;

   assert(scip != NULL);

   exprhdlrdata = SCIPexprhdlrGetData(sourceexprhdlr);
   assert(exprhdlrdata != NULL);
   assert(exprhdlrdata->objexprhdlr != NULL);
   assert(exprhdlrdata->objexprhdlr->scip_ != scip);

   if( exprhdlrdata->objexprhdlr->iscloneable() )
   {
      SCIP_Bool valid = FALSE;
      scip::ObjExprhdlr* newobjexprhdlr;
      newobjexprhdlr = dynamic_cast<scip::ObjExprhdlr*> (exprhdlrdata->objexprhdlr->clone(scip, &valid));

      if( !valid )
         return SCIP_NOMEMORY;

      /* call include method of expression handler object */
      SCIP_EXPRHDLR* targetexprhdlr;
      SCIP_CALL( SCIPincludeObjExprhdlr(scip, &targetexprhdlr, newobjexprhdlr, TRUE) );
   }

   return SCIP_OKAY;
}

/** destructor of expression handler to free user data (called when SCIP is exiting) */
static
SCIP_DECL_EXPRFREEHDLR(exprFreehdlrObj)
{  /*lint --e{715}*/
   assert(exprhdlrdata != NULL);
   assert(*exprhdlrdata != NULL);
   assert((*exprhdlrdata)->objexprhdlr != NULL);
   assert((*exprhdlrdata)->objexprhdlr->scip_ == scip);

   /* call virtual method of exprhdlr object */
   SCIP_CALL( (*exprhdlrdata)->objexprhdlr->scip_freehdlr(scip, exprhdlr,
                                                          exprhdlrdata) );

   /* free exprhdlr object */
   if( (*exprhdlrdata)->deleteobject )
      delete (*exprhdlrdata)->objexprhdlr;

   /* free exprhdlr data */
   delete (*exprhdlrdata);
   *exprhdlrdata = NULL;

   return SCIP_OKAY;
}


/** point evaluation callback of expression handler */
static
SCIP_DECL_EXPREVAL(exprEvalObj)
{  /*lint --e{715}*/
   SCIP_EXPRHDLR* exprhdlr = SCIPexprGetHdlr(expr);
   SCIP_EXPRHDLRDATA* exprhdlrdata = SCIPexprhdlrGetData(exprhdlr);
   assert(exprhdlrdata != NULL);
   assert(exprhdlrdata->objexprhdlr != NULL);

   /* call virtual method of exprhdlr object */
   SCIP_CALL( exprhdlrdata->objexprhdlr->scip_eval(scip, expr, val, sol) );

   return SCIP_OKAY;
}


/** data copy callback of expression handler */
static
SCIP_DECL_EXPRCOPYDATA(exprCopydataObj)
{  /*lint --e{715}*/
   SCIP_EXPRHDLRDATA* exprhdlrdata = SCIPexprhdlrGetData(targetexprhdlr);
   assert(exprhdlrdata != NULL);
   assert(exprhdlrdata->objexprhdlr != NULL);

   /* call virtual method of exprhdlr object */
   SCIP_CALL( exprhdlrdata->objexprhdlr->scip_copydata(targetscip, targetexprhdlr, targetexprdata, sourcescip, sourceexpr) );

   return SCIP_OKAY;
}


/** data free callback of expression handler */
static
SCIP_DECL_EXPRFREEDATA(exprFreedataObj)
{  /*lint --e{715}*/
   SCIP_EXPRHDLR* exprhdlr = SCIPexprGetHdlr(expr);
   SCIP_EXPRHDLRDATA* exprhdlrdata = SCIPexprhdlrGetData(exprhdlr);
   assert(exprhdlrdata != NULL);
   assert(exprhdlrdata->objexprhdlr != NULL);

   /* call virtual method of exprhdlr object */
   SCIP_CALL( exprhdlrdata->objexprhdlr->scip_freedata(scip, expr) );

   return SCIP_OKAY;
}


/** simplify callback of expression handler */
static
SCIP_DECL_EXPRSIMPLIFY(exprSimplifyObj)
{  /*lint --e{715}*/
   SCIP_EXPRHDLR* exprhdlr = SCIPexprGetHdlr(expr);
   SCIP_EXPRHDLRDATA* exprhdlrdata = SCIPexprhdlrGetData(exprhdlr);
   assert(exprhdlrdata != NULL);
   assert(exprhdlrdata->objexprhdlr != NULL);

   /* call virtual method of exprhdlr object */
   SCIP_CALL( exprhdlrdata->objexprhdlr->scip_simplify(scip, expr,
                                                       simplifiedexpr,
                                                       ownercreate,
                                                       ownercreatedata) );

   return SCIP_OKAY;
}


/** compare callback of expression handler */
static
SCIP_DECL_EXPRCOMPARE(exprCompareObj)
{  /*lint --e{715}*/
   SCIP_EXPRHDLR* exprhdlr = SCIPexprGetHdlr(expr1);
   SCIP_EXPRHDLRDATA* exprhdlrdata = SCIPexprhdlrGetData(exprhdlr);
   assert(exprhdlrdata != NULL);
   assert(exprhdlrdata->objexprhdlr != NULL);

   /* call virtual method of exprhdlr object */
   return exprhdlrdata->objexprhdlr->scip_compare(scip, expr1, expr2);
}


/** print callback of expression handler */
static
SCIP_DECL_EXPRPRINT(exprPrintObj)
{  /*lint --e{715}*/
   SCIP_EXPRHDLR* exprhdlr = SCIPexprGetHdlr(expr);
   SCIP_EXPRHDLRDATA* exprhdlrdata = SCIPexprhdlrGetData(exprhdlr);
   assert(exprhdlrdata != NULL);
   assert(exprhdlrdata->objexprhdlr != NULL);

   /* call virtual method of exprhdlr object */
   SCIP_CALL( exprhdlrdata->objexprhdlr->scip_print(scip, expr, stage,
                                                    currentchild,
                                                    parentprecedence, file) );

   return SCIP_OKAY;
}


/** parse callback of expression handler */
static
SCIP_DECL_EXPRPARSE(exprParseObj)
{  /*lint --e{715}*/
   SCIP_EXPRHDLRDATA* exprhdlrdata = SCIPexprhdlrGetData(exprhdlr);
   assert(exprhdlrdata != NULL);
   assert(exprhdlrdata->objexprhdlr != NULL);

   /* call virtual method of exprhdlr object */
   SCIP_CALL( exprhdlrdata->objexprhdlr->scip_parse(scip, exprhdlr, string,
                                                    endstring, expr, success,
                                                    ownercreate,
                                                    ownercreatedata) );

   return SCIP_OKAY;
}


/** backward derivative evaluation callback of expression handler */
static
SCIP_DECL_EXPRBWDIFF(exprBwdiffObj)
{  /*lint --e{715}*/
   SCIP_EXPRHDLR* exprhdlr = SCIPexprGetHdlr(expr);
   SCIP_EXPRHDLRDATA* exprhdlrdata = SCIPexprhdlrGetData(exprhdlr);
   assert(exprhdlrdata != NULL);
   assert(exprhdlrdata->objexprhdlr != NULL);

   /* call virtual method of exprhdlr object */
   SCIP_CALL( exprhdlrdata->objexprhdlr->scip_bwdiff(scip, expr, childidx, val) );

   return SCIP_OKAY;
}


/** forward derivative evaluation callback of expression handler */
static
SCIP_DECL_EXPRFWDIFF(exprFwdiffObj)
{  /*lint --e{715}*/
   SCIP_EXPRHDLR* exprhdlr = SCIPexprGetHdlr(expr);
   SCIP_EXPRHDLRDATA* exprhdlrdata = SCIPexprhdlrGetData(exprhdlr);
   assert(exprhdlrdata != NULL);
   assert(exprhdlrdata->objexprhdlr != NULL);

   /* call virtual method of exprhdlr object */
   SCIP_CALL( exprhdlrdata->objexprhdlr->scip_fwdiff(scip, expr, dot, direction) );

   return SCIP_OKAY;
}


/** backward over forward derivative evaluation callback of expression handler */
static
SCIP_DECL_EXPRBWFWDIFF(exprBwfwdiffObj)
{  /*lint --e{715}*/
   SCIP_EXPRHDLR* exprhdlr = SCIPexprGetHdlr(expr);
   SCIP_EXPRHDLRDATA* exprhdlrdata = SCIPexprhdlrGetData(exprhdlr);
   assert(exprhdlrdata != NULL);
   assert(exprhdlrdata->objexprhdlr != NULL);

   /* call virtual method of exprhdlr object */
   SCIP_CALL( exprhdlrdata->objexprhdlr->scip_bwfwdiff(scip, expr, childidx,
                                                       bardot, direction) );

   return SCIP_OKAY;
}


/** interval evaluation callback of expression handler */
static
SCIP_DECL_EXPRINTEVAL(exprIntevalObj)
{  /*lint --e{715}*/
   SCIP_EXPRHDLR* exprhdlr = SCIPexprGetHdlr(expr);
   SCIP_EXPRHDLRDATA* exprhdlrdata = SCIPexprhdlrGetData(exprhdlr);
   assert(exprhdlrdata != NULL);
   assert(exprhdlrdata->objexprhdlr != NULL);

   /* call virtual method of exprhdlr object */
   SCIP_CALL( exprhdlrdata->objexprhdlr->scip_inteval(scip, expr, interval,
                                                      intevalvar,
                                                      intevalvardata) );

   return SCIP_OKAY;
}


/** estimation callback of expression handler */
static
SCIP_DECL_EXPRESTIMATE(exprEstimateObj)
{  /*lint --e{715}*/
   SCIP_EXPRHDLR* exprhdlr = SCIPexprGetHdlr(expr);
   SCIP_EXPRHDLRDATA* exprhdlrdata = SCIPexprhdlrGetData(exprhdlr);
   assert(exprhdlrdata != NULL);
   assert(exprhdlrdata->objexprhdlr != NULL);

   /* call virtual method of exprhdlr object */
   SCIP_CALL( exprhdlrdata->objexprhdlr->scip_estimate(scip, expr, localbounds,
                                                       globalbounds, refpoint,
                                                       overestimate,
                                                       targetvalue, coefs,
                                                       constant, islocal,
                                                       success, branchcand) );

   return SCIP_OKAY;
}


/** initial estimators callback of expression handler */
static
SCIP_DECL_EXPRINITESTIMATES(exprInitestimatesObj)
{  /*lint --e{715}*/
   SCIP_EXPRHDLR* exprhdlr = SCIPexprGetHdlr(expr);
   SCIP_EXPRHDLRDATA* exprhdlrdata = SCIPexprhdlrGetData(exprhdlr);
   assert(exprhdlrdata != NULL);
   assert(exprhdlrdata->objexprhdlr != NULL);

   /* call virtual method of exprhdlr object */
   SCIP_CALL( exprhdlrdata->objexprhdlr->scip_initestimates(scip, expr, bounds,
                                                            overestimate, coefs,
                                                            constant,
                                                            nreturned) );

   return SCIP_OKAY;
}


/** reverse propagation callback of expression handler */
static
SCIP_DECL_EXPRREVERSEPROP(exprReversepropObj)
{  /*lint --e{715}*/
   SCIP_EXPRHDLR* exprhdlr = SCIPexprGetHdlr(expr);
   SCIP_EXPRHDLRDATA* exprhdlrdata = SCIPexprhdlrGetData(exprhdlr);
   assert(exprhdlrdata != NULL);
   assert(exprhdlrdata->objexprhdlr != NULL);

   /* call virtual method of exprhdlr object */
   SCIP_CALL( exprhdlrdata->objexprhdlr->scip_reverseprop(scip, expr, bounds,
                                                          childrenbounds,
                                                          infeasible) );

   return SCIP_OKAY;
}


/** hash callback of expression handler */
static
SCIP_DECL_EXPRHASH(exprHashObj)
{  /*lint --e{715}*/
   SCIP_EXPRHDLR* exprhdlr = SCIPexprGetHdlr(expr);
   SCIP_EXPRHDLRDATA* exprhdlrdata = SCIPexprhdlrGetData(exprhdlr);
   assert(exprhdlrdata != NULL);
   assert(exprhdlrdata->objexprhdlr != NULL);

   /* call virtual method of exprhdlr object */
   SCIP_CALL( exprhdlrdata->objexprhdlr->scip_hash(scip, expr, hashkey,
                                                   childrenhashes) );

   return SCIP_OKAY;
}


/** curvature callback of expression handler */
static
SCIP_DECL_EXPRCURVATURE(exprCurvatureObj)
{  /*lint --e{715}*/
   SCIP_EXPRHDLR* exprhdlr = SCIPexprGetHdlr(expr);
   SCIP_EXPRHDLRDATA* exprhdlrdata = SCIPexprhdlrGetData(exprhdlr);
   assert(exprhdlrdata != NULL);
   assert(exprhdlrdata->objexprhdlr != NULL);

   /* call virtual method of exprhdlr object */
   SCIP_CALL( exprhdlrdata->objexprhdlr->scip_curvature(scip, expr,
                                                        exprcurvature, success,
                                                        childcurv) );

   return SCIP_OKAY;
}


/** monotonicity callback of expression handler */
static
SCIP_DECL_EXPRMONOTONICITY(exprMonotonicityObj)
{  /*lint --e{715}*/
   SCIP_EXPRHDLR* exprhdlr = SCIPexprGetHdlr(expr);
   SCIP_EXPRHDLRDATA* exprhdlrdata = SCIPexprhdlrGetData(exprhdlr);
   assert(exprhdlrdata != NULL);
   assert(exprhdlrdata->objexprhdlr != NULL);

   /* call virtual method of exprhdlr object */
   SCIP_CALL( exprhdlrdata->objexprhdlr->scip_monotonicity(scip, expr, childidx,
                                                           result) );

   return SCIP_OKAY;
}


/** integrality callback of expression handler */
static
SCIP_DECL_EXPRINTEGRALITY(exprIntegralityObj)
{  /*lint --e{715}*/
   SCIP_EXPRHDLR* exprhdlr = SCIPexprGetHdlr(expr);
   SCIP_EXPRHDLRDATA* exprhdlrdata = SCIPexprhdlrGetData(exprhdlr);
   assert(exprhdlrdata != NULL);
   assert(exprhdlrdata->objexprhdlr != NULL);

   /* call virtual method of exprhdlr object */
   SCIP_CALL( exprhdlrdata->objexprhdlr->scip_integrality(scip, expr,
                                                          integrality) );

   return SCIP_OKAY;
}


/** symmetry information callback of expression handler */
static
SCIP_DECL_EXPRGETSYMDATA(exprGetsymdataObj)
{  /*lint --e{715}*/
   SCIP_EXPRHDLR* exprhdlr = SCIPexprGetHdlr(expr);
   SCIP_EXPRHDLRDATA* exprhdlrdata = SCIPexprhdlrGetData(exprhdlr);
   assert(exprhdlrdata != NULL);
   assert(exprhdlrdata->objexprhdlr != NULL);

   /* call virtual method of exprhdlr object */
   SCIP_CALL( exprhdlrdata->objexprhdlr->scip_getsymdata(scip, expr, symdata) );

   return SCIP_OKAY;
}
}


/*
 * expression handler specific interface methods
 */

/** creates the expression handler for the given expression handler object and includes it in SCIP */
SCIP_RETCODE SCIPincludeObjExprhdlr(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EXPRHDLR**       exprhdlr,           /**< pointer to store the expression handler */   scip::ObjExprhdlr*    objexprhdlr,        /**< expression handler object */
   SCIP_Bool             deleteobject        /**< should the expression handler object be deleted when exprhdlr is freed? */
   )
{
   SCIP_EXPRHDLRDATA* exprhdlrdata;

   assert(scip != NULL);
   assert(objexprhdlr != NULL);
   assert(objexprhdlr->scip_ == scip);

   /* create obj expression handler data */
   exprhdlrdata = new SCIP_EXPRHDLRDATA;
   exprhdlrdata->objexprhdlr = objexprhdlr;
   exprhdlrdata->deleteobject = deleteobject;

   /* include expression handler */
   SCIP_CALL( SCIPincludeExprhdlr(scip, exprhdlr, objexprhdlr->scip_name_,
                                  objexprhdlr->scip_desc_,
                                  objexprhdlr->scip_precedence_, exprEvalObj,
                                  exprhdlrdata) ); /*lint !e429*/
   SCIPexprhdlrSetCopyFreeHdlr(*exprhdlr, exprCopyhdlrObj, exprFreehdlrObj);
   if( objexprhdlr->scip_has_copydata_ || objexprhdlr->scip_has_freedata_ )
      SCIPexprhdlrSetCopyFreeData(*exprhdlr,
                                  objexprhdlr->scip_has_copydata_ ? exprCopydataObj : NULL, objexprhdlr->scip_has_freedata_ ? exprFreedataObj : NULL);
   if( objexprhdlr->scip_has_simplify_ )
      SCIPexprhdlrSetSimplify(*exprhdlr, exprSimplifyObj);
   if( objexprhdlr->scip_has_compare_ )
      SCIPexprhdlrSetCompare(*exprhdlr, exprCompareObj);
   if( objexprhdlr->scip_has_print_ )
      SCIPexprhdlrSetPrint(*exprhdlr, exprPrintObj);
   if( objexprhdlr->scip_has_parse_ )
      SCIPexprhdlrSetParse(*exprhdlr, exprParseObj);
   if( objexprhdlr->scip_has_bwdiff_ || objexprhdlr->scip_has_fwdiff_
      || objexprhdlr->scip_has_bwfwdiff_ )
      SCIPexprhdlrSetDiff(*exprhdlr,
                          objexprhdlr->scip_has_bwdiff_ ? exprBwdiffObj : NULL,
                          objexprhdlr->scip_has_fwdiff_ ? exprFwdiffObj : NULL,
                          objexprhdlr->scip_has_bwfwdiff_ ? exprBwfwdiffObj : NULL);
   if( objexprhdlr->scip_has_inteval_ )
      SCIPexprhdlrSetIntEval(*exprhdlr, exprIntevalObj);
   if( objexprhdlr->scip_has_estimate_ || objexprhdlr->scip_has_initestimates_ )
      SCIPexprhdlrSetEstimate(*exprhdlr,
                              objexprhdlr->scip_has_initestimates_ ? exprInitestimatesObj : NULL,
                              objexprhdlr->scip_has_estimate_ ? exprEstimateObj : NULL);
   if( objexprhdlr->scip_has_reverseprop_ )
      SCIPexprhdlrSetReverseProp(*exprhdlr, exprReversepropObj);
   if( objexprhdlr->scip_has_hash_ )
      SCIPexprhdlrSetHash(*exprhdlr, exprHashObj);
   if( objexprhdlr->scip_has_curvature_ )
      SCIPexprhdlrSetCurvature(*exprhdlr, exprCurvatureObj);
   if( objexprhdlr->scip_has_monotonicity_ )
      SCIPexprhdlrSetMonotonicity(*exprhdlr, exprMonotonicityObj);
   if( objexprhdlr->scip_has_integrality_ )
      SCIPexprhdlrSetIntegrality(*exprhdlr, exprIntegralityObj);
   if( objexprhdlr->scip_has_getsymdata_ )
      SCIPexprhdlrSetGetSymdata(*exprhdlr, exprGetsymdataObj);

   return SCIP_OKAY; /*lint !e429*/
}

/** returns the exprhdlr object of the given name, or 0 if not existing */
scip::ObjExprhdlr* SCIPfindObjExprhdlr(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           name                /**< name of expression handler */
   )
{
   SCIP_EXPRHDLR* exprhdlr;
   SCIP_EXPRHDLRDATA* exprhdlrdata;

   exprhdlr = SCIPfindExprhdlr(scip, name);
   if( exprhdlr == NULL )
      return 0;

   exprhdlrdata = SCIPexprhdlrGetData(exprhdlr);
   assert(exprhdlrdata != NULL);

   return exprhdlrdata->objexprhdlr;
}

/** returns the exprhdlr object for the given expression handler */
scip::ObjExprhdlr* SCIPgetObjExprhdlr(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EXPRHDLR*        exprhdlr            /**< expression handler */
   )
{
   SCIP_EXPRHDLRDATA* exprhdlrdata;

   assert(scip != NULL);
   exprhdlrdata = SCIPexprhdlrGetData(exprhdlr);
   assert(exprhdlrdata != NULL);

   return exprhdlrdata->objexprhdlr;
}
