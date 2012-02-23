/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2012 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License.             */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   objpricer.h
 * @brief  C++ wrapper for variable pricers
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_OBJPRICER_H__
#define __SCIP_OBJPRICER_H__

#include <cstring>

#include "scip/scip.h"
#include "objscip/objprobcloneable.h"

namespace scip
{

/** @brief C++ wrapper for variable pricer
 *
 *  This class defines the interface for variable pricer implemented in C++. Note that there is a pure virtual
 *  function (this function has to be implemented). This function is: scip_redcost().
 *
 *  - \ref PRICER "Instructions for implementing a variable pricer"
 *  - \ref type_pricer.h "Corresponding C interface"
 */
class ObjPricer : public ObjProbCloneable
{
public:
   /*lint --e{1540}*/

   /** SCIP data structure */
   SCIP* scip_;

   /** name of the variable pricer */
   char* scip_name_;
   
   /** description of the variable pricer */
   char* scip_desc_;
   
   /** default priority of the variable pricer */
   const int scip_priority_;

   /** should the pricer be delayed until no other pricers or already existing problem variables with negative reduced
    *  costs are found?
    */
   const SCIP_Bool scip_delay_;

   /** default constructor */
   ObjPricer(
      SCIP*              scip,               /**< SCIP data structure */
      const char*        name,               /**< name of variable pricer */
      const char*        desc,               /**< description of variable pricer */
      int                priority,           /**< priority of the variable pricer */
      SCIP_Bool          delay               /**< should the pricer be delayed until no other pricers or already existing
                                              *   problem variables with negative reduced costs are found?
                                              *   if this is set to FALSE it may happen that the pricer produces columns
                                              *   that already exist in the problem (which are also priced in by the
                                              *   default problem variable pricing in the same round)
                                              */
      )
      : scip_(scip),
        scip_name_(0),
        scip_desc_(0),
        scip_priority_(priority),
        scip_delay_(delay)
   {
      /* the macro SCIPduplicateMemoryArray does not need the first argument: */
      SCIP_CALL_ABORT( SCIPduplicateMemoryArray(scip_, &scip_name_, name, std::strlen(name)+1) );
      SCIP_CALL_ABORT( SCIPduplicateMemoryArray(scip_, &scip_desc_, desc, std::strlen(desc)+1) );
   }

   /** destructor */
   virtual ~ObjPricer()
   {
      /* the macro SCIPfreeMemoryArray does not need the first argument: */
      /*lint --e{64}*/
      SCIPfreeMemoryArray(scip_, &scip_name_);
      SCIPfreeMemoryArray(scip_, &scip_desc_);
   }

   /** destructor of variable pricer to free user data (called when SCIP is exiting) */
   virtual SCIP_RETCODE scip_free(
      SCIP*              scip,               /**< SCIP data structure */
      SCIP_PRICER*       pricer              /**< the variable pricer itself */
      )
   {  /*lint --e{715}*/
      return SCIP_OKAY;
   }
   
   /** initialization method of variable pricer (called after problem was transformed) */
   virtual SCIP_RETCODE scip_init(
      SCIP*              scip,               /**< SCIP data structure */
      SCIP_PRICER*       pricer              /**< the variable pricer itself */
      )
   {  /*lint --e{715}*/
      return SCIP_OKAY;
   }
   
   /** deinitialization method of variable pricer (called before transformed problem is freed) */
   virtual SCIP_RETCODE scip_exit(
      SCIP*              scip,               /**< SCIP data structure */
      SCIP_PRICER*       pricer              /**< the variable pricer itself */
      )
   {  /*lint --e{715}*/
      return SCIP_OKAY;
   }
   
   /** solving process initialization method of variable pricer (called when branch and bound process is about to begin)
    *
    *  This method is called when the presolving was finished and the branch and bound process is about to begin.
    *  The variable pricer may use this call to initialize its branch and bound specific data.
    *
    */
   virtual SCIP_RETCODE scip_initsol(
      SCIP*              scip,               /**< SCIP data structure */
      SCIP_PRICER*       pricer              /**< the variable pricer itself */
      )
   {  /*lint --e{715}*/
      return SCIP_OKAY;
   }
   
   /** solving process deinitialization method of variable pricer (called before branch and bound process data is freed)
    *
    *  This method is called before the branch and bound process is freed.
    *  The variable pricer should use this call to clean up its branch and bound data.
    */
   virtual SCIP_RETCODE scip_exitsol(
      SCIP*              scip,               /**< SCIP data structure */
      SCIP_PRICER*       pricer              /**< the variable pricer itself */
      )
   {  /*lint --e{715}*/
      return SCIP_OKAY;
   }
   
   /** reduced cost pricing method of variable pricer for feasible LPs
    *
    *  Searches for variables that can contribute to improve the current LP's solution value.
    *  In standard branch-and-price, these are variables with negative feasibility, that is negative
    *  reduced costs for non-negative variables, positive reduced costs for non-positive variables,
    *  and non-zero reduced costs for variables that can be negative and positive.
    *
    *  The method is called in the LP solving loop after an LP was proven to be feasible.
    *
    *  Whenever the pricer finds a variable with negative feasibility, it should call SCIPcreateVar()
    *  and SCIPaddPricedVar() to add the variable to the problem. Furthermore, it should call the appropriate
    *  methods of the constraint handlers to add the necessary variable entries to the constraints.
    *
    *  In the usual case that the pricer either adds a new variable or ensures that there are no further variables with negative dual feasibility,
    *  the result pointer should be set to SCIP_SUCCESS. Only if the pricer aborts pricing without creating a new variable, but
    *  there might exist additional variables with negative dual feasibility, the result pointer should be set to SCIP_DIDNOTRUN.
    *  In this case, which sometimes is referred to as "early branching", the lp solution will not be used as a lower bound. 
    *  The pricer can, however, store a valid lower bound in the lowerbound pointer.
    *  If you use your own branching rule (e.g., to branch on constraints), make sure that it is able to branch on pseudo solutions. 
    *  Otherwise, SCIP will use its default branching rules (which all branch on variables). This
    *  could disturb the pricing problem or branching might not even be possible, e.g., if all yet created variables have already been fixed.
    *
    *  input:
    *  - scip            : SCIP main data structure
    *  - pricer          : the variable pricer itself
    *  - lowerbound      : pointer to store a lower bound found by the pricer
    *  - result          : pointer to store the result of the pricer call
    *
    *  possible return values for *result:
    *  - SCIP_SUCCESS    : at least one improving variable was found, or it is ensured that no such variable exists
    *  - SCIP_DIDNOTRUN  : the pricing process was aborted by the pricer, there is no guarantee that the current LP solution is 
    *                      optimal
    *
    */
   virtual SCIP_RETCODE scip_redcost(
      SCIP*              scip,               /**< SCIP data structure */
      SCIP_PRICER*       pricer,             /**< the variable pricer itself */
      SCIP_Real*         lowerbound,         /**< pointer to store a lower bound found by the pricer */
      SCIP_RESULT*       result              /**< pointer to store the result of the pricer call */
      ) = 0;
   
   /** farkas pricing method of variable pricer for infeasible LPs
    *
    *  Searches for variables that can contribute to the feasibility of the current LP.
    *  In standard branch-and-price, these are variables with positive farkas values:
    *
    *  The LP was proven infeasible, so we have an infeasibility proof by the dual farkas multipliers y.
    *  With the values of y, an implicit inequality  y^T A x >= y^T b  is associated, with b given
    *  by the sides of the LP rows and the sign of y:
    *   - if y_i is positive, b_i is the left hand side of the row,
    *   - if y_i is negative, b_i is the right hand side of the row.
    *
    *  y is chosen in a way, such that the valid inequality  y^T A x >= y^T b  is violated by all x,
    *  especially by the (for this inequality least infeasible solution) x' defined by 
    *     x'_i := ub_i, if y^T A_i >= 0
    *     x'_i := lb_i, if y^T A_i < 0.
    *  Pricing in this case means to add variables i with positive farkas value, i.e. y^T A_i x'_i > 0.
    *
    *  The method is called in the LP solving loop after an LP was proven to be infeasible.
    *
    *  Whenever the pricer finds a variable with positive farkas value, it should call SCIPcreateVar()
    *  and SCIPaddPricedVar() to add the variable to the problem. Furthermore, it should call the appropriate
    *  methods of the constraint handlers to add the necessary variable entries to the constraints.
    */
   virtual SCIP_RETCODE scip_farkas(
      SCIP*              scip,               /**< SCIP data structure */
      SCIP_PRICER*       pricer              /**< the variable pricer itself */
      )
   {  /*lint --e{715}*/
      return SCIP_OKAY;
   }
};

} /* namespace scip */


   
/** creates the variable pricer for the given variable pricer object and includes it in SCIP
 *
 *  The method should be called in one of the following ways:
 *
 *   1. The user is resposible of deleting the object:
 *       SCIP_CALL( SCIPcreate(&scip) );
 *       ...
 *       MyPricer* mypricer = new MyPricer(...);
 *       SCIP_CALL( SCIPincludeObjPricer(scip, &mypricer, FALSE) );
 *       ...
 *       SCIP_CALL( SCIPfree(&scip) );
 *       delete mypricer;    // delete pricer AFTER SCIPfree() !
 *
 *   2. The object pointer is passed to SCIP and deleted by SCIP in the SCIPfree() call:
 *       SCIP_CALL( SCIPcreate(&scip) );
 *       ...
 *       SCIP_CALL( SCIPincludeObjPricer(scip, new MyPricer(...), TRUE) );
 *       ...
 *       SCIP_CALL( SCIPfree(&scip) );  // destructor of MyPricer is called here
 */
extern
SCIP_RETCODE SCIPincludeObjPricer(
   SCIP*                 scip,               /**< SCIP data structure */
   scip::ObjPricer*      objpricer,          /**< variable pricer object */
   SCIP_Bool             deleteobject        /**< should the pricer object be deleted when pricer is freed? */
   );

/** returns the variable pricer object of the given name, or 0 if not existing */
extern
scip::ObjPricer* SCIPfindObjPricer(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           name                /**< name of variable pricer */
   );
   
/** returns the variable pricer object for the given pricer */
extern
scip::ObjPricer* SCIPgetObjPricer(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PRICER*          pricer              /**< pricer */
   );

#endif
