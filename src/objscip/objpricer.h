/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2004 Tobias Achterberg                              */
/*                                                                           */
/*                  2002-2004 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the SCIP Academic Licence.        */
/*                                                                           */
/*  You should have received a copy of the SCIP Academic License             */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#pragma ident "@(#) $Id: objpricer.h,v 1.9 2004/12/14 12:08:01 bzfpfend Exp $"

/**@file   objpricer.h
 * @brief  C++ wrapper for variable pricers
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __OBJPRICER_H__
#define __OBJPRICER_H__


extern "C" 
{
#include "scip.h"
}


namespace scip
{

/** C++ wrapper object for variable pricers */
class ObjPricer
{
public:
   /** name of the variable pricer */
   const char* const scip_name_;
   
   /** description of the variable pricer */
   const char* const scip_desc_;
   
   /** default priority of the variable pricer */
   const int scip_priority_;

   /** default constructor */
   ObjPricer(
      const char*   name,               /**< name of variable pricer */
      const char*   desc,               /**< description of variable pricer */
      int           priority            /**< priority of the variable pricer */
      )
      : scip_name_(name),
        scip_desc_(desc),
        scip_priority_(priority)
   {
   }

   /** destructor */
   virtual ~ObjPricer()
   {
   }

   /** destructor of variable pricer to free user data (called when SCIP is exiting) */
   virtual RETCODE scip_free(
      SCIP*         scip,               /**< SCIP data structure */
      PRICER*       pricer              /**< the variable pricer itself */
      )
   {
      return SCIP_OKAY;
   }
   
   /** initialization method of variable pricer (called after problem was transformed) */
   virtual RETCODE scip_init(
      SCIP*         scip,               /**< SCIP data structure */
      PRICER*       pricer              /**< the variable pricer itself */
      )
   {
      return SCIP_OKAY;
   }
   
   /** deinitialization method of variable pricer (called before transformed problem is freed) */
   virtual RETCODE scip_exit(
      SCIP*         scip,               /**< SCIP data structure */
      PRICER*       pricer              /**< the variable pricer itself */
      )
   {
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
    */
   virtual RETCODE scip_redcost(
      SCIP*         scip,               /**< SCIP data structure */
      PRICER*       pricer              /**< the variable pricer itself */
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
   virtual RETCODE scip_farkas(
      SCIP*         scip,               /**< SCIP data structure */
      PRICER*       pricer              /**< the variable pricer itself */
      )
   {
      return SCIP_OKAY;
   }
};

} /* namespace scip */


   
/** creates the variable pricer for the given variable pricer object and includes it in SCIP
 *
 *  The method should be called in one of the following ways:
 *
 *   1. The user is resposible of deleting the object:
 *       CHECK_OKAY( SCIPcreate(&scip) );
 *       ...
 *       MyPricer* mypricer = new MyPricer(...);
 *       CHECK_OKAY( SCIPincludeObjPricer(scip, &mypricer, FALSE) );
 *       ...
 *       CHECK_OKAY( SCIPfree(&scip) );
 *       delete mypricer;    // delete pricer AFTER SCIPfree() !
 *
 *   2. The object pointer is passed to SCIP and deleted by SCIP in the SCIPfree() call:
 *       CHECK_OKAY( SCIPcreate(&scip) );
 *       ...
 *       CHECK_OKAY( SCIPincludeObjPricer(scip, new MyPricer(...), TRUE) );
 *       ...
 *       CHECK_OKAY( SCIPfree(&scip) );  // destructor of MyPricer is called here
 */
extern
RETCODE SCIPincludeObjPricer(
   SCIP*            scip,               /**< SCIP data structure */
   scip::ObjPricer* objpricer,          /**< variable pricer object */
   Bool             deleteobject        /**< should the pricer object be deleted when pricer is freed? */
   );

/** returns the variable pricer object of the given name, or NULL if not existing */
extern
scip::ObjPricer* SCIPfindObjPricer(
   SCIP*            scip,               /**< SCIP data structure */
   const char*      name                /**< name of variable pricer */
   );
   
/** returns the variable pricer object for the given pricer */
extern
scip::ObjPricer* SCIPgetObjPricer(
   SCIP*            scip,               /**< SCIP data structure */
   PRICER*          pricer              /**< pricer */
   );

#endif
