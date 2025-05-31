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

/**@file   objexprhdlr.h
 * @brief  C++ wrapper for expression handlers
 * @author Kevin Kofler
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_OBJEXPRHDLR_H__
#define __SCIP_OBJEXPRHDLR_H__


#include <cassert>
#include <cstring>
#include <utility>

#include "scip/scip.h"
#include "objscip/objprobcloneable.h"

namespace scip
{

/** @brief C++ wrapper for expression handlers
 *
 *  This class defines the interface for expression handlers implemented in C++. Note that there is a pure virtual
 *  function (which has to be implemented): the function scip_eval().
 *
 *  - \ref EXPRHDLR "Instructions for implementing an expression handler"
 *  - \ref EXPRHDLRS "List of available expression handlers"
 *  - \ref type_expr.h "Corresponding C interface"
 */
class ObjExprhdlr : public ObjProbCloneable
{
public:
   /*lint --e{1540}*/

   /** SCIP data structure */
   SCIP* scip_;

   /** name of the expression handler */
   char* scip_name_;

   /** description of the expression handler */
   char* scip_desc_;

   /** precedence of expression operation relative to other expression (used for printing) */
   const unsigned int scip_precedence_;

   /** whether scip_copydata is implemented */
   const SCIP_Bool scip_has_copydata_;

   /** whether scip_freedata is implemented */
   const SCIP_Bool scip_has_freedata_;

   /** whether scip_simplify is implemented */
   const SCIP_Bool scip_has_simplify_;

   /** whether scip_compare is implemented */
   const SCIP_Bool scip_has_compare_;

   /** whether scip_print is implemented */
   const SCIP_Bool scip_has_print_;

   /** whether scip_parse is implemented */
   const SCIP_Bool scip_has_parse_;

   /** whether scip_bwdiff is implemented */
   const SCIP_Bool scip_has_bwdiff_;

   /** whether scip_fwdiff is implemented */
   const SCIP_Bool scip_has_fwdiff_;

   /** whether scip_bwfwdiff is implemented */
   const SCIP_Bool scip_has_bwfwdiff_;

   /** whether scip_inteval is implemented */
   const SCIP_Bool scip_has_inteval_;

   /** whether scip_estimate is implemented */
   const SCIP_Bool scip_has_estimate_;

   /** whether scip_initestimates is implemented */
   const SCIP_Bool scip_has_initestimates_;

   /** whether scip_reverseprop is implemented */
   const SCIP_Bool scip_has_reverseprop_;

   /** whether scip_hash is implemented */
   const SCIP_Bool scip_has_hash_;

   /** whether scip_curvature is implemented */
   const SCIP_Bool scip_has_curvature_;

   /** whether scip_monotonicity is implemented */
   const SCIP_Bool scip_has_monotonicity_;

   /** whether scip_integrality is implemented */
   const SCIP_Bool scip_has_integrality_;

   /** whether scip_getsymdata is implemented */
   const SCIP_Bool scip_has_getsymdata_;

   /** default constructor */
   ObjExprhdlr(
      SCIP*              scip,              /**< SCIP data structure */
      const char*        name,              /**< name of expression handler */
      const char*        desc,              /**< description of expression handler */
      unsigned int       precedence,        /**< precedence of expression operation */
      SCIP_Bool          has_copydata,      /**< whether scip_copydata is implemented */
      SCIP_Bool          has_freedata,      /**< whether scip_freedata is implemented */
      SCIP_Bool          has_simplify,      /**< whether scip_simplify is implemented */
      SCIP_Bool          has_compare,       /**< whether scip_compare is implemented */
      SCIP_Bool          has_print,         /**< whether scip_print is implemented */
      SCIP_Bool          has_parse,         /**< whether scip_parse is implemented */
      SCIP_Bool          has_bwdiff,        /**< whether scip_bwdiff is implemented */
      SCIP_Bool          has_fwdiff,        /**< whether scip_fwdiff is implemented */
      SCIP_Bool          has_bwfwdiff,      /**< whether scip_bwfwdiff is implemented */
      SCIP_Bool          has_inteval,       /**< whether scip_inteval is implemented */
      SCIP_Bool          has_estimate,      /**< whether scip_estimate is implemented */
      SCIP_Bool          has_initestimates, /**< whether scip_initestimates is implemented */
      SCIP_Bool          has_reverseprop,   /**< whether scip_reverseprop is implemented */
      SCIP_Bool          has_hash,          /**< whether scip_hash is implemented */
      SCIP_Bool          has_curvature,     /**< whether scip_curvature is implemented */
      SCIP_Bool          has_monotonicity,  /**< whether scip_monotonicity is implemented */
      SCIP_Bool          has_integrality,   /**< whether scip_integrality is implemented */
      SCIP_Bool          has_getsymdata     /**< whether scip_getsymdata is implemented */
      )
      : scip_(scip),
        scip_name_(0),
        scip_desc_(0),
        scip_precedence_(precedence),
        scip_has_copydata_(has_copydata),
        scip_has_freedata_(has_freedata),
        scip_has_simplify_(has_simplify),
        scip_has_compare_(has_compare),
        scip_has_print_(has_print),
        scip_has_parse_(has_parse),
        scip_has_bwdiff_(has_bwdiff),
        scip_has_fwdiff_(has_fwdiff),
        scip_has_bwfwdiff_(has_bwfwdiff),
        scip_has_inteval_(has_inteval),
        scip_has_estimate_(has_estimate),
        scip_has_initestimates_(has_initestimates),
        scip_has_reverseprop_(has_reverseprop),
        scip_has_hash_(has_hash),
        scip_has_curvature_(has_curvature),
        scip_has_monotonicity_(has_monotonicity),
        scip_has_integrality_(has_integrality),
        scip_has_getsymdata_(has_getsymdata)
   {
      /* the macro SCIPduplicateMemoryArray does not need the first argument: */
      SCIP_CALL_ABORT( SCIPduplicateMemoryArray(scip_, &scip_name_, name, std::strlen(name)+1) );
      SCIP_CALL_ABORT( SCIPduplicateMemoryArray(scip_, &scip_desc_, desc, std::strlen(desc)+1) );
   }

   /** copy constructor */
   ObjExprhdlr(const ObjExprhdlr& o)
       : ObjExprhdlr(o.scip_, o.scip_name_, o.scip_desc_, o.scip_precedence_,
                     o.scip_has_copydata_, o.scip_has_freedata_,
                     o.scip_has_simplify_, o.scip_has_compare_,
                     o.scip_has_print_, o.scip_has_parse_, o.scip_has_bwdiff_,
                     o.scip_has_fwdiff_, o.scip_has_bwfwdiff_,
                     o.scip_has_inteval_, o.scip_has_estimate_,
                     o.scip_has_initestimates_, o.scip_has_reverseprop_,
                     o.scip_has_hash_, o.scip_has_curvature_,
                     o.scip_has_monotonicity_, o.scip_has_integrality_,
                     o.scip_has_getsymdata_)
   {
   }

   /** move constructor */
   ObjExprhdlr(ObjExprhdlr&& o)
       : scip_(o.scip_),
         scip_name_(0),
         scip_desc_(0),
         scip_precedence_(o.scip_precedence_),
         scip_has_copydata_(o.scip_has_copydata_),
         scip_has_freedata_(o.scip_has_freedata_),
         scip_has_simplify_(o.scip_has_simplify_),
         scip_has_compare_(o.scip_has_compare_),
         scip_has_print_(o.scip_has_print_),
         scip_has_parse_(o.scip_has_parse_),
         scip_has_bwdiff_(o.scip_has_bwdiff_),
         scip_has_fwdiff_(o.scip_has_fwdiff_),
         scip_has_bwfwdiff_(o.scip_has_bwfwdiff_),
         scip_has_inteval_(o.scip_has_inteval_),
         scip_has_estimate_(o.scip_has_estimate_),
         scip_has_initestimates_(o.scip_has_initestimates_),
         scip_has_reverseprop_(o.scip_has_reverseprop_),
         scip_has_hash_(o.scip_has_hash_),
         scip_has_curvature_(o.scip_has_curvature_),
         scip_has_monotonicity_(o.scip_has_monotonicity_),
         scip_has_integrality_(o.scip_has_integrality_),
         scip_has_getsymdata_(o.scip_has_getsymdata_)
   {
      std::swap(scip_name_, o.scip_name_);
      std::swap(scip_desc_, o.scip_desc_);
   }

   /** destructor */
   virtual ~ObjExprhdlr()
   {
      /* the macro SCIPfreeMemoryArray does not need the first argument: */
      /*lint --e{64}*/
      SCIPfreeMemoryArray(scip_, &scip_name_);
      SCIPfreeMemoryArray(scip_, &scip_desc_);
   }

   /** assignment of polymorphic classes causes slicing and is therefore disabled. */
   ObjExprhdlr& operator=(const ObjExprhdlr& o) = delete;

   /** assignment of polymorphic classes causes slicing and is therefore disabled. */
   ObjExprhdlr& operator=(ObjExprhdlr&& o) = delete;

   /** destructor of expression handler to free user data (called when SCIP is exiting)
    *
    *  @see SCIP_DECL_EXPRFREEHDLR(x) in @ref type_expr.h
    */
   virtual SCIP_DECL_EXPRFREEHDLR(scip_freehdlr)
   {  /*lint --e{715}*/
      return SCIP_OKAY;
   }

   /** point evaluation callback of expression handler
    *
    *  @see SCIP_DECL_EXPREVAL(x) in @ref type_expr.h
    */
   virtual SCIP_DECL_EXPREVAL(scip_eval) = 0;

   /** data copy callback of expression handler
    *
    *  This method MUST be overridden if scip_has_copydata_ is TRUE.
    *
    *  @see SCIP_DECL_EXPRCOPYDATA(x) in @ref type_expr.h
    */
   virtual SCIP_DECL_EXPRCOPYDATA(scip_copydata)
   {  /*lint --e{715}*/
      /* This method MUST be overridden if scip_has_copydata_ is TRUE. */
      return SCIP_NOTIMPLEMENTED;
   }

   /** data free callback of expression handler
    *
    *  This method MUST be overridden if scip_has_freedata_ is TRUE.
    *
    *  @see SCIP_DECL_EXPRFREEDATA(x) in @ref type_expr.h
    */
   virtual SCIP_DECL_EXPRFREEDATA(scip_freedata)
   {  /*lint --e{715}*/
      /* This method MUST be overridden if scip_has_freedata_ is TRUE. */
      return SCIP_NOTIMPLEMENTED;
   }

   /** simplify callback of expression handler
    *
    *  This method MUST be overridden if scip_has_simplify_ is TRUE.
    *
    *  @see SCIP_DECL_EXPRSIMPLIFY(x) in @ref type_expr.h
    */
   virtual SCIP_DECL_EXPRSIMPLIFY(scip_simplify)
   {  /*lint --e{715}*/
      /* This method MUST be overridden if scip_has_simplify_ is TRUE. */
      return SCIP_NOTIMPLEMENTED;
   }

   /** compare callback of expression handler
    *
    *  This method MUST be overridden if scip_has_compare_ is TRUE.
    *
    *  @see SCIP_DECL_EXPRCOMPARE(x) in @ref type_expr.h
    */
   virtual SCIP_DECL_EXPRCOMPARE(scip_compare)
   {  /*lint --e{715}*/
      /* This method MUST be overridden if scip_has_compare_ is TRUE. */
      return 0;
   }

   /** print callback of expression handler
    *
    *  This method MUST be overridden if scip_has_print_ is TRUE.
    *
    *  @see SCIP_DECL_EXPRPRINT(x) in @ref type_expr.h
    */
   virtual SCIP_DECL_EXPRPRINT(scip_print)
   {  /*lint --e{715}*/
      /* This method MUST be overridden if scip_has_print_ is TRUE. */
      return SCIP_NOTIMPLEMENTED;
   }

   /** parse callback of expression handler
    *
    *  This method MUST be overridden if scip_has_parse_ is TRUE.
    *
    *  @see SCIP_DECL_EXPRPARSE(x) in @ref type_expr.h
    */
   virtual SCIP_DECL_EXPRPARSE(scip_parse)
   {  /*lint --e{715}*/
      /* This method MUST be overridden if scip_has_parse_ is TRUE. */
      return SCIP_NOTIMPLEMENTED;
   }

   /** backward derivative evaluation callback of expression handler
    *
    *  This method MUST be overridden if scip_has_bwdiff_ is TRUE.
    *
    *  @see SCIP_DECL_EXPRBWDIFF(x) in @ref type_expr.h
    */
   virtual SCIP_DECL_EXPRBWDIFF(scip_bwdiff)
   {  /*lint --e{715}*/
      /* This method MUST be overridden if scip_has_bwdiff_ is TRUE. */
      return SCIP_NOTIMPLEMENTED;
   }

   /** forward derivative evaluation callback of expression handler
    *
    *  This method MUST be overridden if scip_has_fwdiff_ is TRUE.
    *
    *  @see SCIP_DECL_EXPRFWDIFF(x) in @ref type_expr.h
    */
   virtual SCIP_DECL_EXPRFWDIFF(scip_fwdiff)
   {  /*lint --e{715}*/
      /* This method MUST be overridden if scip_has_fwdiff_ is TRUE. */
      return SCIP_NOTIMPLEMENTED;
   }

   /** backward over forward derivative evaluation callback of expression handler
    *
    *  This method MUST be overridden if scip_has_bwfwdiff_ is TRUE.
    *
    *  @see SCIP_DECL_EXPRBWFWDIFF(x) in @ref type_expr.h
    */
   virtual SCIP_DECL_EXPRBWFWDIFF(scip_bwfwdiff)
   {  /*lint --e{715}*/
      /* This method MUST be overridden if scip_has_bwfwdiff_ is TRUE. */
      return SCIP_NOTIMPLEMENTED;
   }

   /** interval evaluation callback of expression handler
    *
    *  This method MUST be overridden if scip_has_inteval_ is TRUE.
    *
    *  @see SCIP_DECL_EXPRINTEVAL(x) in @ref type_expr.h
    */
   virtual SCIP_DECL_EXPRINTEVAL(scip_inteval)
   {  /*lint --e{715}*/
      /* This method MUST be overridden if scip_has_inteval_ is TRUE. */
      return SCIP_NOTIMPLEMENTED;
   }

   /** estimation callback of expression handler
    *
    *  This method MUST be overridden if scip_has_estimate_ is TRUE.
    *
    *  @see SCIP_DECL_EXPRESTIMATE(x) in @ref type_expr.h
    */
   virtual SCIP_DECL_EXPRESTIMATE(scip_estimate)
   {  /*lint --e{715}*/
      /* This method MUST be overridden if scip_has_estimate_ is TRUE. */
      return SCIP_NOTIMPLEMENTED;
   }

   /** initial estimators callback of expression handler
    *
    *  This method MUST be overridden if scip_has_initestimates_ is TRUE.
    *
    *  @see SCIP_DECL_EXPRINITESTIMATES(x) in @ref type_expr.h
    */
   virtual SCIP_DECL_EXPRINITESTIMATES(scip_initestimates)
   {  /*lint --e{715}*/
      /* This method MUST be overridden if scip_has_initestimates_ is TRUE. */
      return SCIP_NOTIMPLEMENTED;
   }

   /** reverse propagation callback of expression handler
    *
    *  This method MUST be overridden if scip_has_reverseprop_ is TRUE.
    *
    *  @see SCIP_DECL_EXPRREVERSEPROP(x) in @ref type_expr.h
    */
   virtual SCIP_DECL_EXPRREVERSEPROP(scip_reverseprop)
   {  /*lint --e{715}*/
      /* This method MUST be overridden if scip_has_reverseprop_ is TRUE. */
      return SCIP_NOTIMPLEMENTED;
   }

   /** hash callback of expression handler
    *
    *  This method MUST be overridden if scip_has_hash_ is TRUE.
    *
    *  @see SCIP_DECL_EXPRHASH(x) in @ref type_expr.h
    */
   virtual SCIP_DECL_EXPRHASH(scip_hash)
   {  /*lint --e{715}*/
      /* This method MUST be overridden if scip_has_hash_ is TRUE. */
      return SCIP_NOTIMPLEMENTED;
   }

   /** curvature callback of expression handler
    *
    *  This method MUST be overridden if scip_has_curvature_ is TRUE.
    *
    *  @see SCIP_DECL_EXPRCURVATURE(x) in @ref type_expr.h
    */
   virtual SCIP_DECL_EXPRCURVATURE(scip_curvature)
   {  /*lint --e{715}*/
      /* This method MUST be overridden if scip_has_curvature_ is TRUE. */
      return SCIP_NOTIMPLEMENTED;
   }

   /** monotonicity callback of expression handler
    *
    *  This method MUST be overridden if scip_has_monotonicity_ is TRUE.
    *
    *  @see SCIP_DECL_EXPRMONOTONICITY(x) in @ref type_expr.h
    */
   virtual SCIP_DECL_EXPRMONOTONICITY(scip_monotonicity)
   {  /*lint --e{715}*/
      /* This method MUST be overridden if scip_has_monotonicity_ is TRUE. */
      return SCIP_NOTIMPLEMENTED;
   }

   /** integrality callback of expression handler
    *
    *  This method MUST be overridden if scip_has_integrality_ is TRUE.
    *
    *  @see SCIP_DECL_EXPRINTEGRALITY(x) in @ref type_expr.h
    */
   virtual SCIP_DECL_EXPRINTEGRALITY(scip_integrality)
   {  /*lint --e{715}*/
      /* This method MUST be overridden if scip_has_integrality_ is TRUE. */
      return SCIP_NOTIMPLEMENTED;
   }

   /** symmetry information callback of expression handler
    *
    *  This method MUST be overridden if scip_has_getsymdata_ is TRUE.
    *
    *  @see SCIP_DECL_EXPRGETSYMDATA(x) in @ref type_expr.h
    */
   virtual SCIP_DECL_EXPRGETSYMDATA(scip_getsymdata)
   {  /*lint --e{715}*/
      /* This method MUST be overridden if scip_has_getsymdata_ is TRUE. */
      return SCIP_NOTIMPLEMENTED;
   }
};

} /* namespace scip */



/** creates the expression handler for the given expression handler object and includes it in SCIP
 *
 *  The method should be called in one of the following ways:
 *
 *   1. The user is resposible of deleting the object:
 *       SCIP_CALL( SCIPcreate(&scip) );
 *       ...
 *       SCIP_EXPRHDLR* cexprhdlr;
 *       MyExprhdlr* myexprhdlr = new MyExprhdlr(...);
 *       SCIP_CALL( SCIPincludeObjExprhdlr(scip, &cexprhdlr, &myexprhdlr, FALSE) );
 *       ...
 *       SCIP_CALL( SCIPfree(&scip) );
 *       delete myexprhdlr;    // delete exprhdlr AFTER SCIPfree() !
 *
 *   2. The object pointer is passed to SCIP and deleted by SCIP in the SCIPfree() call:
 *       SCIP_CALL( SCIPcreate(&scip) );
 *       ...
 *       SCIP_EXPRHDLR* cexprhdlr;
 *       SCIP_CALL( SCIPincludeObjExprhdlr(scip, &cexprhdlr, new MyExprhdlr(...), TRUE) );
 *       ...
 *       SCIP_CALL( SCIPfree(&scip) );  // destructor of MyExprhdlr is called here
 */
SCIP_EXPORT
SCIP_RETCODE SCIPincludeObjExprhdlr(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EXPRHDLR**       exprhdlr,           /**< pointer to store the expression handler */   scip::ObjExprhdlr*    objconshdlr,        /**< expression handler object */
   SCIP_Bool             deleteobject        /**< should the expression handler object be deleted when exprhdlr is freed? */
   );

/** returns the exprhdlr object of the given name, or 0 if not existing */
SCIP_EXPORT
scip::ObjExprhdlr* SCIPfindObjExprhdlr(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           name                /**< name of expression handler */
   );

/** returns the exprhdlr object for the given expression handler */
SCIP_EXPORT
scip::ObjExprhdlr* SCIPgetObjExprhdlr(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EXPRHDLR*        exprhdlr            /**< expression handler */
   );

#endif
