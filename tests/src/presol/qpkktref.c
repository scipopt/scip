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

/**@file   unittest-qpkktref.c
 * @brief  Unittest for adding the KKT conditions to a quadratic program
 * @author Tobias Fischer
 */

/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>

#include "scip/scip.h"
#include "scip/scipshell.h"
#include "scip/scipdefplugins.h"
#include "scip/dialog_default.h"



/** for a linear constraint \f$a^T x \leq b\f$, create the complementarity constraint \f$\mu \cdot s = 0\f$, where
 *  \f$s = b - a^T x\f$ and \f$\mu\f$ is the dual variable associated to the constraint \f$a^T x \leq b\f$
 */
static
SCIP_RETCODE createKKTComplementarityLinear(
   SCIP*                 scip,               /**< SCIP pointer */
   const char*           namepart,           /**< name of linear constraint */
   SCIP_VAR**            vars,               /**< variables of linear constraint */
   SCIP_Real*            vals,               /**< coefficients of variables in linear constraint */
   SCIP_Real             lhs,                /**< left hand side of linear constraint */
   SCIP_Real             rhs,                /**< right hand side of linear constraint */
   int                   nvars,              /**< number of variables of linear constraint */
   SCIP_VAR*             dualvar,            /**< dual variable associated to linear constraint */
   SCIP_Bool             takelhs,            /**< whether to consider the lhs or the rhs of the constraint */
   int*                  naddconss           /**< buffer to increase with number of created additional constraints */
   )
{
   char name[SCIP_MAXSTRLEN];
   SCIP_CONS* KKTlincons;
   SCIP_CONS* sos1cons;
   SCIP_VAR* slack;
   SCIP_Real slackcoef;
   SCIP_Real eqval;

   assert( scip != NULL );
   assert( namepart != NULL );
   assert( vars != NULL );
   assert( vals != NULL );
   assert( dualvar != NULL );
   assert( ! takelhs || ! SCIPisInfinity(scip, -lhs) );
   assert( takelhs || ! SCIPisInfinity(scip, rhs) );
   assert( naddconss != NULL );

   if ( takelhs )
   {
      eqval = lhs;
      slackcoef = -1.0;
      (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "slack_lhs_%s", namepart);
   }
   else
   {
      eqval = rhs;
      slackcoef = 1.0;
      (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "slack_rhs_%s", namepart);
   }

   /* create slack variable */
   SCIP_CALL( SCIPcreateVarBasic(scip, &slack, name, 0.0, SCIPinfinity(scip), 0.0, SCIP_VARTYPE_CONTINUOUS) );

   /* add skack variable */
   SCIP_CALL( SCIPaddVar(scip, slack) );

   /* create a new linear constraint */
   (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "KKTlin_%s_%d", namepart, takelhs);
   SCIP_CALL( SCIPcreateConsBasicLinear(scip, &KKTlincons, name, nvars, vars, vals, eqval, eqval) );

   /* add slack variable to linear constraint */
   SCIP_CALL( SCIPaddCoefLinear(scip, KKTlincons, slack, slackcoef) );

   /* create SOS1 (complementarity) constraint involving dual variable of linear constraint and slack variable */
   (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "KKTsos1_lin_%s_%d", namepart, takelhs);
   SCIP_CALL( SCIPcreateConsBasicSOS1(scip, &sos1cons, name, 0, NULL, NULL) );

   /* add slack and dual variable to SOS1 constraint */
   SCIP_CALL( SCIPaddVarSOS1(scip, sos1cons, slack, 1.0) );
   SCIP_CALL( SCIPaddVarSOS1(scip, sos1cons, dualvar, 2.0) );

   /* add/release constraints */
   SCIP_CALL( SCIPaddCons(scip, sos1cons) );
   SCIP_CALL( SCIPaddCons(scip, KKTlincons) );
   SCIP_CALL( SCIPreleaseCons(scip, &sos1cons) );
   SCIP_CALL( SCIPreleaseCons(scip, &KKTlincons) );
   *naddconss = *naddconss + 2;

   /* release slack variable */
   SCIP_CALL( SCIPreleaseVar(scip, &slack) );

   return SCIP_OKAY;
}

/** create complementarity constraints of KKT conditions associated to bounds of variables
 * - for an upper bound constraint \f$x_i \leq u_i\f$, create the complementarity constraint \f$\mu_i \cdot s_i = 0\f$,
 *   where \f$s_i = u_i - x_i\f$ and \f$\mu_i\f$ is the dual variable of the upper bound constraint
 * - for a lower bound constraint \f$x_i \geq l_i\f$, create the complementarity constraint \f$\lambda_i \cdot w_i = 0\f$,
 *   where \f$w_i = x_i - l_i\f$
 *   and \f$\lambda_i\f$ is the dual variable of the lower bound constraint
 */
static
SCIP_RETCODE createKKTComplementarityBounds(
   SCIP*                 scip,               /**< SCIP pointer */
   SCIP_VAR*             var,                /**< variable */
   SCIP_VAR*             dualvar,            /**< dual variable associated to bound of variable */
   SCIP_Bool             takelb,             /**< whether to consider the lower or upper bound of variable */
   int*                  naddconss           /**< buffer to increase with number of created additional constraints */
   )
{
   char name[SCIP_MAXSTRLEN];
   SCIP_CONS* KKTlincons;
   SCIP_CONS* sos1cons;
   SCIP_VAR* slack;
   SCIP_Real slackcoef;
   SCIP_Real eqval;

   assert( scip != NULL );
   assert( var != NULL );
   assert( dualvar != NULL );
   assert( ! takelb || ! SCIPisInfinity(scip, -SCIPvarGetLbGlobal(var)) );
   assert( takelb || ! SCIPisInfinity(scip, SCIPvarGetUbGlobal(var)) );
   assert( naddconss != NULL );

   if ( takelb )
   {
      eqval = SCIPvarGetLbGlobal(var);
      slackcoef = -1.0;
      (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "slack_lb_%s", SCIPvarGetName(var));
   }
   else
   {
      eqval = SCIPvarGetUbGlobal(var);
      slackcoef = 1.0;
      (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "slack_ub_%s", SCIPvarGetName(var));
   }

   /* create complementarity constraint; if bound is nonzero, we additionally need to introduce a slack variable */
   if ( SCIPisFeasZero(scip, eqval) && SCIPvarGetStatus(var) != SCIP_VARSTATUS_MULTAGGR )
   {
      /* create SOS1 (complementarity) constraint involving dual variable of linear constraint and slack variable */
      (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "KKTsos1_bound%s_%d", SCIPvarGetName(var), takelb);
      SCIP_CALL( SCIPcreateConsBasicSOS1(scip, &sos1cons, name, 0, NULL, NULL) );

      /* add slack and dual variable to SOS1 constraint */
      SCIP_CALL( SCIPaddVarSOS1(scip, sos1cons, var, 1.0) );
      SCIP_CALL( SCIPaddVarSOS1(scip, sos1cons, dualvar, 2.0) );

      /* add/release constraint */
      SCIP_CALL( SCIPaddCons(scip, sos1cons) );
      SCIP_CALL( SCIPreleaseCons(scip, &sos1cons) );
      ++(*naddconss);
   }
   else
   {
      /* create slack variable */
      SCIP_CALL( SCIPcreateVarBasic(scip, &slack, name, 0.0, SCIPinfinity(scip), 0.0, SCIP_VARTYPE_CONTINUOUS) );

      /* add skack variable */
      SCIP_CALL( SCIPaddVar(scip, slack) );

      /* create a new linear constraint */
      (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "KKT_bound%s_%d", SCIPvarGetName(var), takelb);
      SCIP_CALL( SCIPcreateConsBasicLinear(scip, &KKTlincons, name, 0, NULL, NULL, eqval, eqval) );

      /* add slack variable to linear constraint */
      SCIP_CALL( SCIPaddCoefLinear(scip, KKTlincons, var, 1.0) );
      SCIP_CALL( SCIPaddCoefLinear(scip, KKTlincons, slack, slackcoef) );

      /* create SOS1 (complementarity) constraint involving dual variable of linear constraint and slack variable */
      (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "KKTsos1_bound%s_%d", SCIPvarGetName(var), takelb);
      SCIP_CALL( SCIPcreateConsBasicSOS1(scip, &sos1cons, name, 0, NULL, NULL) );

      /* add slack and dual variable to SOS1 constraint */
      SCIP_CALL( SCIPaddVarSOS1(scip, sos1cons, slack, 1.0) );
      SCIP_CALL( SCIPaddVarSOS1(scip, sos1cons, dualvar, 2.0) );

      /* add/release constraints */
      SCIP_CALL( SCIPaddCons(scip, sos1cons) );
      SCIP_CALL( SCIPaddCons(scip, KKTlincons) );
      SCIP_CALL( SCIPreleaseCons(scip, &sos1cons) );
      SCIP_CALL( SCIPreleaseCons(scip, &KKTlincons) );
      *naddconss = *naddconss + 2;

      /* release slack variable */
      SCIP_CALL( SCIPreleaseVar(scip, &slack) );
   }

   return SCIP_OKAY;
}

/** create the complementarity constraints of the KKT-like conditions associated to a binary variable \f$x_i\f$;
 *  these are \f$(1 - x_i) \cdot z_i = 0\f$ and \f$x_i \cdot (z_i - \lambda_i) = 0\f$, where \f$z_i\f$ and
 *  \f$\lambda_i\f$ are dual variables
 */
static
SCIP_RETCODE createKKTComplementarityBinary(
   SCIP*                 scip,               /**< SCIP pointer */
   SCIP_VAR*             var,                /**< variable */
   SCIP_VAR*             dualbin1,           /**< first dual variable associated to binary variable */
   SCIP_VAR*             dualbin2,           /**< second dual variable associated to binary variable */
   int*                  naddconss           /**< buffer to increase with number of created additional constraints */
   )
{
   char name[SCIP_MAXSTRLEN];
   SCIP_CONS* conslinbin1;
   SCIP_CONS* conslinbin2;
   SCIP_CONS* sos1cons1;
   SCIP_CONS* sos1cons2;
   SCIP_VAR* slackbin1;
   SCIP_VAR* slackbin2;

   assert( scip != NULL );
   assert( var != NULL );
   assert( dualbin1 != NULL );
   assert( dualbin2 != NULL );
   assert( naddconss != NULL );

   /* create first slack variable associated to binary constraint; domain [-inf, inf] */
   (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "dual_%s_slackbin1", SCIPvarGetName(var));
   SCIP_CALL( SCIPcreateVarBasic(scip, &slackbin1, name, -SCIPinfinity(scip), SCIPinfinity(scip), 0.0,
        SCIP_VARTYPE_CONTINUOUS) );
   SCIP_CALL( SCIPaddVar(scip, slackbin1) );
   assert( slackbin1 != NULL );

   /* create a new linear constraint: dualbin1 - dualbin2 = slackbin */
   (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "KKTBinary1_%s", SCIPvarGetName(var));
   SCIP_CALL( SCIPcreateConsBasicLinear(scip, &conslinbin1, name, 0, NULL, NULL, 0.0, 0.0) );
   SCIP_CALL( SCIPaddCoefLinear(scip, conslinbin1, dualbin1, 1.0) );
   SCIP_CALL( SCIPaddCoefLinear(scip, conslinbin1, dualbin2, -1.0) );
   SCIP_CALL( SCIPaddCoefLinear(scip, conslinbin1, slackbin1, -1.0) );
   SCIP_CALL( SCIPaddCons(scip, conslinbin1) );
   SCIP_CALL( SCIPreleaseCons(scip, &conslinbin1) );
   ++(*naddconss);

   /* create SOS1 (complementarity) constraint involving binary variable and slack variable */
   (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "KKTsos1_bin1%s", SCIPvarGetName(var));
   SCIP_CALL( SCIPcreateConsBasicSOS1(scip, &sos1cons1, name, 0, NULL, NULL) );

   /* add slack and dual variable to SOS1 constraint */
   SCIP_CALL( SCIPaddVarSOS1(scip, sos1cons1, var, 1.0) );
   SCIP_CALL( SCIPaddVarSOS1(scip, sos1cons1, slackbin1, 2.0) );

   /* add/release constraint */
   SCIP_CALL( SCIPaddCons(scip, sos1cons1) );
   SCIP_CALL( SCIPreleaseCons(scip, &sos1cons1) );
   ++(*naddconss);

   /* release slack variable */
   SCIP_CALL( SCIPreleaseVar(scip, &slackbin1) );


   /* create second slack variable associated to binary constraint; domain [0, inf] */
   (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "dual_%s_slackbin2", SCIPvarGetName(var));
   SCIP_CALL( SCIPcreateVarBasic(scip, &slackbin2, name, 0.0, SCIPinfinity(scip), 0.0, SCIP_VARTYPE_CONTINUOUS) );
   SCIP_CALL( SCIPaddVar(scip, slackbin2) );
   assert( slackbin2 != NULL );

   /* create a new linear constraint: 1.0 - var = slackbin2 */
   (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "KKTBinary2_%s", SCIPvarGetName(var));
   SCIP_CALL( SCIPcreateConsBasicLinear(scip, &conslinbin2, name, 0, NULL, NULL, 1.0, 1.0) );
   SCIP_CALL( SCIPaddCoefLinear(scip, conslinbin2, var, 1.0) );
   SCIP_CALL( SCIPaddCoefLinear(scip, conslinbin2, slackbin2, 1.0) );
   SCIP_CALL( SCIPaddCons(scip, conslinbin2) );
   SCIP_CALL( SCIPreleaseCons(scip, &conslinbin2) );
   ++(*naddconss);

   /* create SOS1 (complementarity) constraint involving first dual variable and slack variable */
   (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "KKTsos1_bin2%s", SCIPvarGetName(var));
   SCIP_CALL( SCIPcreateConsBasicSOS1(scip, &sos1cons2, name, 0, NULL, NULL) );

   /* add slack and dual variable to SOS1 constraint */
   SCIP_CALL( SCIPaddVarSOS1(scip, sos1cons2, dualbin1, 1.0) );
   SCIP_CALL( SCIPaddVarSOS1(scip, sos1cons2, slackbin2, 2.0) );

   /* add/release constraint */
   SCIP_CALL( SCIPaddCons(scip, sos1cons2) );
   SCIP_CALL( SCIPreleaseCons(scip, &sos1cons2) );
   ++(*naddconss);

   /* release slack variable */
   SCIP_CALL( SCIPreleaseVar(scip, &slackbin2) );

   return SCIP_OKAY;
}

/** create/get dual constraint of KKT conditions associated to primal variable @n@n
 * if variable does not already exist in hashmap then
 * 1. create dual constraint for variable
 * 2. create a dual variable \f$\mu_i\f$ for the upper bound constraint \f$x_i \leq u_i\f$
 * 3. create a dual variable \f$\lambda_i\f$ for the lower bound constraint \f$x_i \geq l_i\f$
 * 4. create the complementarity constraint \f$\mu_i \cdot s_i = 0\f$, where \f$s_i = u_i - x_i\f$
 * 5. create the complementarity constraint \f$\lambda_i \cdot w_i = 0\f$, where \f$w_i = x_i - l_i\f$
 * 6. add objective coefficients of dual variables
 * 7. the treatment of binary variables needs special care see the documentation of createKKTComplementarityBinary()
 *
 * if variable exists in hasmap then the dual constraint associated to the variable has already been created and is returned
 */
static
SCIP_RETCODE createKKTDualCons(
   SCIP*                 scip,               /**< SCIP pointer */
   SCIP_CONS*            objcons,            /**< objective constraint */
   SCIP_VAR*             var,                /**< variable */
   SCIP_HASHMAP*         varhash,            /**< hash map from variable to index of linear constraint */
   SCIP_CONS**           dualconss,          /**< array with dual constraints */
   int*                  ndualconss,         /**< pointer to store number of dual constraints */
   SCIP_CONS**           dualcons,           /**< dual constraint associated to variable */
   int*                  naddconss           /**< buffer to increase with number of created additional constraints */
   )
{
   SCIP_VAR* dualub = NULL;     /* dual variable associated to upper bound constraint */
   SCIP_VAR* duallb = NULL;     /* dual variable associated to lower bound constraint */
   SCIP_VAR* dualbin1 = NULL;   /* first dual variable associated to binary variable */
   SCIP_VAR* dualbin2 = NULL;   /* second dual variable associated to binary variable */

   assert( scip != NULL );
   assert( objcons != NULL );
   assert( var != NULL );
   assert( varhash != NULL );
   assert( dualconss != NULL );
   assert( ndualconss != NULL );
   assert( naddconss != NULL );

   /* if variable exists in hashmap */
   if ( SCIPhashmapExists(varhash, var) )
   {
      int ind;
      ind = (int) (size_t) SCIPhashmapGetImage(varhash, var);
      *dualcons = dualconss[ind];
   }
   else
   {
      char name[SCIP_MAXSTRLEN];
      SCIP_Real lb;
      SCIP_Real ub;

      lb = SCIPvarGetLbGlobal(var);
      ub = SCIPvarGetUbGlobal(var);

      /* create dual variables corresponding to the bounds of the variables; binary variables have to be treated in a
       * different way */
      if ( SCIPvarIsBinary(var) )
      {
         /* create first dual variable associated to binary constraint; the domain of dualbin is [-inf,inf]; the objective
          * coefficient is -0.5 */
         (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "dual_%s_bin1", SCIPvarGetName(var));
         SCIP_CALL( SCIPcreateVarBasic(scip, &dualbin1, name, -SCIPinfinity(scip), SCIPinfinity(scip), 0.0,
              SCIP_VARTYPE_CONTINUOUS) );
         SCIP_CALL( SCIPaddVar(scip, dualbin1) );
         assert( dualbin1 != NULL );
         SCIP_CALL( SCIPaddCoefLinear(scip, objcons, dualbin1, -0.5) );

         /* create second variable associated to binary constraint; the domain of dualbin2 is [-inf,inf]; the objective
          * coefficient is zero */
         (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "dual_%s_bin2", SCIPvarGetName(var));
         SCIP_CALL( SCIPcreateVarBasic(scip, &dualbin2, name, -SCIPinfinity(scip), SCIPinfinity(scip), 0.0,
              SCIP_VARTYPE_CONTINUOUS) );
         SCIP_CALL( SCIPaddVar(scip, dualbin2) );
         assert( dualbin2 != NULL );
      }
      else
      {
         if ( ! SCIPisInfinity(scip, -lb) )
         {
            /* create dual variable associated to lower bound; the domain of duallb is [0,inf] */
            (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "dual_%s_lb", SCIPvarGetName(var));
            SCIP_CALL( SCIPcreateVarBasic(scip, &duallb, name, 0.0, SCIPinfinity(scip), 0.0, SCIP_VARTYPE_CONTINUOUS) );
            SCIP_CALL( SCIPaddVar(scip, duallb) );
            assert( duallb != NULL );
            SCIP_CALL( SCIPaddCoefLinear(scip, objcons, duallb, 0.5 * lb) );
         }

         if ( ! SCIPisInfinity(scip, ub) )
         {
            /* create dual variable associated to upper bound; the domain of dualub is [0,inf] */
            (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "dual_%s_ub", SCIPvarGetName(var));
            SCIP_CALL( SCIPcreateVarBasic(scip, &dualub, name, 0.0, SCIPinfinity(scip), 0.0, SCIP_VARTYPE_CONTINUOUS) );
            SCIP_CALL( SCIPaddVar(scip, dualub) );
            assert( dualub != NULL );
            SCIP_CALL( SCIPaddCoefLinear(scip, objcons, dualub, -0.5 * ub) );
         }
      }

      /* add variable in map  */
      SCIP_CALL( SCIPhashmapInsert(varhash, var, (void*) (size_t) *ndualconss) );/*lint !e571*/
      assert( *ndualconss == (int) (size_t) SCIPhashmapGetImage(varhash, var) );

      /* create a new linear constraint */
      (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "KKTref_%s", SCIPvarGetName(var));
      SCIP_CALL( SCIPcreateConsBasicLinear(scip, dualcons, name, 0, NULL, NULL, 0.0, 0.0) );

      /* add dual constraint to array for later use */
      dualconss[(*ndualconss)++] = *dualcons;


      /* add dual variables to dual constraints and create complementarity constraints; binary variables have to be
       * treated in a different way */
      if ( SCIPvarIsBinary(var) )
      {
         /* add coefficient of second dual variable corresponding to binary variable */
         SCIP_CALL( SCIPaddCoefLinear(scip, *dualcons, dualbin2, 1.0) );

         /* create complementarity constraints */
         SCIP_CALL( createKKTComplementarityBinary(scip, var, dualbin1, dualbin2, naddconss) );

         SCIP_CALL( SCIPreleaseVar(scip, &dualbin1) );
         SCIP_CALL( SCIPreleaseVar(scip, &dualbin2) );
      }
      else
      {
         if ( duallb != NULL )
         {
            /* add dual variable corresponding to lower bound of variable */
            SCIP_CALL( SCIPaddCoefLinear(scip, *dualcons, duallb, -1.0) );

            /* create complementarity constraint between slack variable of lower bound constraint and dual variable of
             * lower bound */
            SCIP_CALL( createKKTComplementarityBounds(scip, var, duallb, TRUE, naddconss) );

            SCIP_CALL( SCIPreleaseVar(scip, &duallb) );
         }

         if ( dualub != NULL )
         {
            /* add dual variable corresponding to upper bound of variable */
            SCIP_CALL( SCIPaddCoefLinear(scip, *dualcons, dualub, 1.0) );

            /* create complementarity constraint between slack variable of upper bound constraint and dual variable of
             * upper bound */
            SCIP_CALL( createKKTComplementarityBounds(scip, var, dualub, FALSE, naddconss) );

            SCIP_CALL( SCIPreleaseVar(scip, &dualub) );
         }
      }
   }
   assert( *dualcons != NULL );

   return SCIP_OKAY;
}

/** handle (a single) linear constraint for quadratic constraint update
 * 1. create the dual constraints (i.e., the two rows of \f$Q x + c + A^T \mu = 0\f$) associated to the variables of the
 *    linear constraint, if not done already
 * 2. create the dual variables and the complementarity constraints for the lower and upper bound constraints of the
 *    variables of the linear constraint, if not done already
 * 3. create the dual variable \f$\mu_i\f$ associated to this linear constraint
 * 4. create the complementarity constraint \f$\mu_i \cdot (Ax - b)_i = 0\f$ associated to this linear constraint
 * 5. add objective coefficients of dual variables
 *
 * for steps 1 and 2 see the documentation of createKKTDualCons() for further information.@n
 * for step 4 see the documentation of the function createKKTComplementarityLinear() for further information.
 */
static
SCIP_RETCODE presolveAddKKTLinearCons(
   SCIP*                 scip,               /**< SCIP pointer */
   SCIP_CONS*            objcons,            /**< objective constraint */
   const char*           namepart,           /**< name of linear constraint */
   SCIP_VAR**            vars,               /**< variables of linear constraint */
   SCIP_Real*            vals,               /**< coefficients of variables in linear constraint */
   SCIP_Real             lhs,                /**< left hand side of linear constraint */
   SCIP_Real             rhs,                /**< right hand side of linear constraint */
   int                   nvars,              /**< number of variables of linear constraint */
   SCIP_HASHMAP*         varhash,            /**< hash map from variable to index of linear constraint */
   SCIP_CONS**           dualconss,          /**< array with dual constraints */
   int*                  ndualconss,         /**< pointer to store number of dual constraints */
   int*                  naddconss           /**< buffer to increase with number of created additional constraints */
   )
{
   int i;

   assert( scip != NULL );
   assert( objcons != NULL );
   assert( namepart != NULL );
   assert( varhash != NULL );
   assert( dualconss != NULL );
   assert( ndualconss != NULL );
   assert( vars != NULL );
   assert( vals != NULL );
   assert( namepart != NULL );
   assert( naddconss != NULL );

   /* differ between left hand side and right hand side case (i=0 -> lhs; i=1 -> rhs) */
   for (i = 0; i < 2; ++i)
   {
      char name[SCIP_MAXSTRLEN];
      SCIP_VAR* duallin = NULL;
      int j;

      /* scip one iteration if lhs equals rhs */
      if ( SCIPisFeasEQ(scip, lhs, rhs) )
         i = 1;

      /* create dual variable corresponding to linear constraint */
      if ( i == 0 )
      {
         assert( ! SCIPisFeasEQ(scip, lhs, rhs) );

         if ( SCIPisInfinity(scip, -lhs) )
            continue;

         (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "dual_%s_lhs", namepart);
         SCIP_CALL( SCIPcreateVarBasic(scip, &duallin, name, 0.0, SCIPinfinity(scip), 0.0, SCIP_VARTYPE_CONTINUOUS) );
         SCIP_CALL( SCIPaddVar(scip, duallin) );
         SCIP_CALL( SCIPaddCoefLinear(scip, objcons, duallin, 0.5 * lhs) );

         /* create complementarity constraint between dual variable and slack variable of linear constraint */
         SCIP_CALL( createKKTComplementarityLinear(scip, namepart, vars, vals, lhs, rhs, nvars, duallin, TRUE, naddconss) );
      }
      else
      {
         if ( SCIPisInfinity(scip, rhs) )
            continue;

         (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "dual_%s_rhs", namepart);
         if ( SCIPisFeasEQ(scip, lhs, rhs) )
         {
            SCIP_CALL( SCIPcreateVarBasic(scip, &duallin, name, -SCIPinfinity(scip), SCIPinfinity(scip), 0.0,
                 SCIP_VARTYPE_CONTINUOUS) );
            SCIP_CALL( SCIPaddVar(scip, duallin) );
            SCIP_CALL( SCIPaddCoefLinear(scip, objcons, duallin, -0.5 * rhs) );
         }
         else
         {
            SCIP_CALL( SCIPcreateVarBasic(scip, &duallin, name, 0.0, SCIPinfinity(scip), 0.0, SCIP_VARTYPE_CONTINUOUS) );
            SCIP_CALL( SCIPaddVar(scip, duallin) );
            SCIP_CALL( SCIPaddCoefLinear(scip, objcons, duallin, -0.5 * rhs) );

            /* create complementarity constraint between dual variable and slack variable of linear constraint */
            SCIP_CALL( createKKTComplementarityLinear(scip, namepart, vars, vals, lhs, rhs, nvars, duallin, FALSE,
                 naddconss) );
         }
      }
      assert( duallin != NULL );


      /* loop through variables of linear constraint */
      for (j = 0; j < nvars; ++j)
      {
         SCIP_CONS* dualcons = NULL;  /* dual constraint associated to variable */
         SCIP_VAR* var;

         var = vars[j];

         /* create/get dual constraint associated to variable;
          * if variable does not already exist in hashmap then create dual variables for its bounds */
         SCIP_CALL( createKKTDualCons(scip, objcons, var, varhash, dualconss, ndualconss, &dualcons, naddconss) );
         assert( dualcons != NULL );

         /* add dual variable corresponding to linear constraint */
         if ( i == 0 )
         {
            SCIP_CALL( SCIPaddCoefLinear(scip, dualcons, duallin, -vals[j]) );
         }
         else
         {
            SCIP_CALL( SCIPaddCoefLinear(scip, dualcons, duallin, vals[j]) );
         }
      }

      /* release dual variable */
      SCIP_CALL( SCIPreleaseVar(scip, &duallin) );
   }

   return SCIP_OKAY;
}


/* check number of added KKT constraints */
static
SCIP_RETCODE checkNAddConss(
   SCIP*                 scip,               /**< SCIP instance */
   SCIP_VAR**            vars,               /**< variables */
   int                   nvars,              /**< number of variables */
   SCIP_Real             lhs,                /**< left hand side of linear constraint */
   SCIP_Real             rhs,                /**< right hand side of linear constraint */
   int                   naddconss           /**< number of added constraints */
   )
{
   int nconss = 0;
   int v;

   for (v = 0; v < nvars; ++v)
   {
      SCIP_VAR* var = vars[v];

      if ( SCIPvarIsBinary(var) )
         nconss += 4;   /* two linear constraints (for slack variables) and two SOS1 constraints */
      else
      {
         SCIP_Real ub;
         SCIP_Real lb;

         lb = SCIPvarGetLbGlobal(var);
         ub = SCIPvarGetUbGlobal(var);

         if ( ! SCIPisFeasZero(scip, lb) && ! SCIPisInfinity(scip, -lb) )
            nconss += 2;   /* one linear (for slack variable) and one SOS1 constraint */

         if ( ! SCIPisFeasZero(scip, ub) && ! SCIPisInfinity(scip, ub) )
            nconss += 2;   /* one linear (for slack variable) and one SOS1 constraint */
      }
   }

   if ( ! SCIPisFeasEQ(scip, lhs, rhs) )
   {
      if ( ! SCIPisInfinity(scip, -lhs) )
      {
         nconss += 2; /* one linear (for slack variable) and one SOS1 constraint */
      }

      if ( ! SCIPisInfinity(scip, rhs) )
      {
         nconss += 2; /* one linear (for slack variable) and one SOS1 constraint */
      }
   }

   if ( naddconss != nconss )
   {
      SCIPerrorMessage("failed test for number of added SOS1 constraints and constraints for slack variables.\n");
      return SCIP_ERROR;
   }

   return SCIP_OKAY;
}


/* check number of dual constraints */
static
SCIP_RETCODE checkNDualConss(
   SCIP*                 scip,               /**< SCIP instance */
   int                   nvars,              /**< number of variables */
   int                   ndualconss          /**< number of dual constraints */
   )
{
   if ( ndualconss != nvars )
   {
      SCIPerrorMessage("failed test for number of dual constraints.\n");
      return SCIP_ERROR;
   }

   return SCIP_OKAY;
}


/** checks linear constraint */
static
SCIP_RETCODE checkConsLinear(
   SCIP*                 scip,               /**< SCIP instance */
   SCIP_CONS*            cons                /**< linear constraint */
   )
{
   SCIP_HASHMAP* varhash; /* hash map from variable to index of dual constraint */
   SCIP_CONS* objcons;
   SCIP_CONS** dualconss;
   SCIP_VAR** vars;
   SCIP_Real* vals;
   SCIP_Real lhs;
   SCIP_Real rhs;
   int ndualconss = 0;
   int naddconss = 0;
   int nvars;
   int j;

   SCIP_CALL( SCIPcreateConsBasicLinear(scip, &objcons, "objcons", 0, NULL, NULL, 0.0, 0.0) );

   /* set up hash map */
   SCIP_CALL( SCIPhashmapCreate(&varhash, SCIPblkmem(scip),
         SCIPcalcHashtableSize(10 * (SCIPgetNVars(scip) + SCIPgetNFixedVars(scip)))) );

   /* allocate buffer array */
   SCIP_CALL( SCIPallocBufferArray(scip, &dualconss, 2 * SCIPgetNVars(scip) + 2 * SCIPgetNFixedVars(scip)) ); /*lint !e647*/

   lhs = SCIPgetLhsLinear(scip, cons);
   rhs = SCIPgetRhsLinear(scip, cons);
   nvars = SCIPgetNVarsLinear(scip, cons);
   vars = SCIPgetVarsLinear(scip, cons);
   vals = SCIPgetValsLinear(scip, cons);

   /* handle linear constraint for quadratic constraint update */
   SCIP_CALL( presolveAddKKTLinearCons(scip, objcons, SCIPconsGetName(cons),
        vars, vals, lhs, rhs, nvars, varhash, dualconss, &ndualconss, &naddconss) );

   /* check number of added constraints */
   SCIP_CALL( checkNDualConss(scip, nvars, ndualconss) );
   SCIP_CALL( checkNAddConss(scip, vars, nvars, lhs, rhs, naddconss) );

   /* free buffer array */
   for (j = 0; j < ndualconss; ++j)
   {
      SCIP_CALL( SCIPreleaseCons(scip, &dualconss[j]) );
   }
   SCIPfreeBufferArray(scip, &dualconss);

   /* free hash map */
   SCIPhashmapFree(&varhash);

   SCIP_CALL( SCIPreleaseCons(scip, &objcons) );

   return SCIP_OKAY;
}


/** checks knapsack constraint */
static
SCIP_RETCODE checkConsKnapsack(
   SCIP*                 scip,               /**< SCIP instance */
   SCIP_CONS*            cons                /**< knapsack constraint */
   )
{
   SCIP_HASHMAP* varhash; /* hash map from variable to index of dual constraint */
   SCIP_CONS* objcons;
   SCIP_CONS** dualconss;
   SCIP_VAR** vars;
   SCIP_Longint* weights;
   SCIP_Real* vals;
   SCIP_Real lhs;
   SCIP_Real rhs;
   int ndualconss = 0;
   int naddconss = 0;
   int nvars;
   int v;
   int j;

   SCIP_CALL( SCIPcreateConsBasicLinear(scip, &objcons, "objcons", 0, NULL, NULL, 0.0, 0.0) );

   /* set up hash map */
   SCIP_CALL( SCIPhashmapCreate(&varhash, SCIPblkmem(scip),
                                SCIPcalcHashtableSize(10 * (SCIPgetNVars(scip) + SCIPgetNFixedVars(scip)))) );

   /* allocate buffer array */
   SCIP_CALL( SCIPallocBufferArray(scip, &dualconss, 2 * SCIPgetNVars(scip) + 2 * SCIPgetNFixedVars(scip)) ); /*lint !e647*/

   lhs = -SCIPinfinity(scip);
   rhs = (SCIP_Real) SCIPgetCapacityKnapsack(scip, cons);
   nvars = SCIPgetNVarsKnapsack(scip, cons);
   vars = SCIPgetVarsKnapsack(scip, cons);
   weights = SCIPgetWeightsKnapsack(scip, cons);

   /* set coefficients of variables */
   SCIP_CALL( SCIPallocBufferArray(scip, &vals, nvars) );
   for (v = 0; v < nvars; ++v)
      vals[v] = (SCIP_Real) weights[v];

   /* handle linear constraint for quadratic constraint update */
   SCIP_CALL( presolveAddKKTLinearCons(scip, objcons, SCIPconsGetName(cons),
         vars, vals, lhs, rhs, nvars, varhash, dualconss, &ndualconss, &naddconss) );

   /* check number of added constraints */
   SCIP_CALL( checkNDualConss(scip, nvars, ndualconss) );
   SCIP_CALL( checkNAddConss(scip, vars, nvars, lhs, rhs, naddconss) );

   /* free buffer arrays */
   for (j = 0; j < ndualconss; ++j)
   {
      SCIP_CALL( SCIPreleaseCons(scip, &dualconss[j]) );
   }
   SCIPfreeBufferArray(scip, &vals);
   SCIPfreeBufferArray(scip, &dualconss);

   /* free hash map */
   SCIPhashmapFree(&varhash);

   SCIP_CALL( SCIPreleaseCons(scip, &objcons) );

   return SCIP_OKAY;
}


/** create quadratic problem with one linear constraint, set objective sense to minimize */
static
SCIP_RETCODE createProbLinear1(
   SCIP*                 scip                /**< SCIP instance */
   )
{
   SCIP_CONS* cons;
   SCIP_VAR* xvar;
   SCIP_VAR* yvar;
   SCIP_VAR* zvar;
   SCIP_VAR* objvar;
   SCIP_VAR* vars[3];
   SCIP_Real vals[3];

   SCIP_CALL( SCIPcreateProbBasic(scip, "qpkktref-linear1") );

   /* create variables */
   SCIP_CALL( SCIPcreateVarBasic(scip, &xvar, "x", 0.0, 1.0, 0.0, SCIP_VARTYPE_BINARY) );
   SCIP_CALL( SCIPcreateVarBasic(scip, &yvar, "y", 0.0, 1.0, 0.0, SCIP_VARTYPE_BINARY) );
   SCIP_CALL( SCIPcreateVarBasic(scip, &zvar, "z", -1.0, 1.0, 0.0, SCIP_VARTYPE_CONTINUOUS) );
   SCIP_CALL( SCIPcreateVarBasic(scip, &objvar, "objvar", -SCIPinfinity(scip), SCIPinfinity(scip), 1.0, SCIP_VARTYPE_CONTINUOUS) );

   SCIP_CALL( SCIPaddVar(scip, xvar) );
   SCIP_CALL( SCIPaddVar(scip, yvar) );
   SCIP_CALL( SCIPaddVar(scip, zvar) );
   SCIP_CALL( SCIPaddVar(scip, objvar) );

   /* create quadratic inequality */
   SCIP_CALL( SCIPcreateConsBasicQuadratic(scip, &cons, "quadratic", 0, NULL, NULL, 0, NULL, NULL, NULL, 2.0, 2.0) );

   SCIP_CALL( SCIPaddLinearVarQuadratic(scip, cons, objvar, 1.0) );
   SCIP_CALL( SCIPaddLinearVarQuadratic(scip, cons, zvar, -3.0) );

   SCIP_CALL( SCIPaddSquareCoefQuadratic(scip, cons, xvar, 2.0) );
   SCIP_CALL( SCIPaddSquareCoefQuadratic(scip, cons, yvar, -1.0) );
   SCIP_CALL( SCIPaddSquareCoefQuadratic(scip, cons, zvar, 1.0) );

   SCIP_CALL( SCIPaddBilinTermQuadratic(scip, cons, xvar, yvar, -2.0) );
   SCIP_CALL( SCIPaddBilinTermQuadratic(scip, cons, xvar, zvar, 1.0) );
   SCIP_CALL( SCIPaddBilinTermQuadratic(scip, cons, zvar, yvar, 1.0) );

   SCIP_CALL( SCIPaddCons(scip, cons) );
   SCIP_CALL( SCIPreleaseCons(scip, &cons) );

   /* create linear inequality */
   vars[0] = xvar;
   vars[1] = yvar;
   vars[2] = zvar;

   vals[0] = 1.0;
   vals[1] = 1.0;
   vals[2] = 1.0;

   SCIP_CALL( SCIPcreateConsBasicLinear(scip, &cons, "lower", 3, vars, vals, 0.25, 0.75) );
   SCIP_CALL( SCIPaddCons(scip, cons) );

   /* check constraint */
   SCIP_CALL( checkConsLinear(scip, cons) );

   SCIP_CALL( SCIPreleaseCons(scip, &cons) );
   SCIP_CALL( SCIPreleaseVar(scip, &xvar) );
   SCIP_CALL( SCIPreleaseVar(scip, &yvar) );
   SCIP_CALL( SCIPreleaseVar(scip, &zvar) );
   SCIP_CALL( SCIPreleaseVar(scip, &objvar) );

   SCIP_CALL( SCIPprintOrigProblem(scip, NULL, NULL, FALSE) );

   return SCIP_OKAY;
}


/** create quadratic problem with one linear constraint, set objective sense to maximize */
static
SCIP_RETCODE createProbLinear2(
   SCIP*                 scip                /**< SCIP instance */
   )
{
   SCIP_CONS* cons;
   SCIP_VAR* xvar;
   SCIP_VAR* yvar;
   SCIP_VAR* zvar;
   SCIP_VAR* objvar;
   SCIP_VAR* vars[3];
   SCIP_Real vals[3];

   SCIP_CALL( SCIPcreateProbBasic(scip, "qpkktref-linear2") );

   SCIP_CALL( SCIPsetObjsense(scip, SCIP_OBJSENSE_MAXIMIZE) );

   /* create variables */
   SCIP_CALL( SCIPcreateVarBasic(scip, &xvar, "x", 0.0, 1.0, 0.0, SCIP_VARTYPE_BINARY) );
   SCIP_CALL( SCIPcreateVarBasic(scip, &yvar, "y", 0.0, 1.0, 0.0, SCIP_VARTYPE_BINARY) );
   SCIP_CALL( SCIPcreateVarBasic(scip, &zvar, "z", -1.0, 1.0, 0.0, SCIP_VARTYPE_CONTINUOUS) );
   SCIP_CALL( SCIPcreateVarBasic(scip, &objvar, "objvar", -SCIPinfinity(scip), SCIPinfinity(scip), 1.0, SCIP_VARTYPE_CONTINUOUS) );

   SCIP_CALL( SCIPaddVar(scip, xvar) );
   SCIP_CALL( SCIPaddVar(scip, yvar) );
   SCIP_CALL( SCIPaddVar(scip, zvar) );
   SCIP_CALL( SCIPaddVar(scip, objvar) );

   /* create quadratic inequality */
   SCIP_CALL( SCIPcreateConsBasicQuadratic(scip, &cons, "quadratic", 0, NULL, NULL, 0, NULL, NULL, NULL, 2.0, 2.0) );

   SCIP_CALL( SCIPaddLinearVarQuadratic(scip, cons, objvar, 1.0) );
   SCIP_CALL( SCIPaddLinearVarQuadratic(scip, cons, zvar, -3.0) );

   SCIP_CALL( SCIPaddSquareCoefQuadratic(scip, cons, xvar, 2.0) );
   SCIP_CALL( SCIPaddSquareCoefQuadratic(scip, cons, yvar, -1.0) );
   SCIP_CALL( SCIPaddSquareCoefQuadratic(scip, cons, zvar, 1.0) );

   SCIP_CALL( SCIPaddBilinTermQuadratic(scip, cons, xvar, yvar, -2.0) );
   SCIP_CALL( SCIPaddBilinTermQuadratic(scip, cons, xvar, zvar, 1.0) );
   SCIP_CALL( SCIPaddBilinTermQuadratic(scip, cons, zvar, yvar, 1.0) );

   SCIP_CALL( SCIPaddCons(scip, cons) );
   SCIP_CALL( SCIPreleaseCons(scip, &cons) );

   /* create linear inequality */
   vars[0] = xvar;
   vars[1] = yvar;
   vars[2] = zvar;

   vals[0] = 1.0;
   vals[1] = 1.0;
   vals[2] = 1.0;

   SCIP_CALL( SCIPcreateConsBasicLinear(scip, &cons, "lower", 3, vars, vals, 0.25, 0.75) );
   SCIP_CALL( SCIPaddCons(scip, cons) );

   /* check constraint */
   SCIP_CALL( checkConsLinear(scip, cons) );

   SCIP_CALL( SCIPreleaseCons(scip, &cons) );
   SCIP_CALL( SCIPreleaseVar(scip, &xvar) );
   SCIP_CALL( SCIPreleaseVar(scip, &yvar) );
   SCIP_CALL( SCIPreleaseVar(scip, &zvar) );
   SCIP_CALL( SCIPreleaseVar(scip, &objvar) );

   SCIP_CALL( SCIPprintOrigProblem(scip, NULL, NULL, FALSE) );

   return SCIP_OKAY;
}


/** create quadratic problem with one knapsack constraint */
static
SCIP_RETCODE createProbKnapsack(
   SCIP*                 scip                /**< SCIP instance */
   )
{
   SCIP_CONS* cons;
   SCIP_VAR* xvar;
   SCIP_VAR* yvar;
   SCIP_VAR* zvar;
   SCIP_VAR* objvar;
   SCIP_VAR* vars[2];
   SCIP_Longint vals[2];

   SCIP_CALL( SCIPcreateProbBasic(scip, "qpkktref-knapsack") );

   SCIP_CALL( SCIPsetObjsense(scip, SCIP_OBJSENSE_MAXIMIZE) );

   /* create variables */
   SCIP_CALL( SCIPcreateVarBasic(scip, &xvar, "x", 0.0, 1.0, 0.0, SCIP_VARTYPE_BINARY) );
   SCIP_CALL( SCIPcreateVarBasic(scip, &yvar, "y", 0.0, 1.0, 0.0, SCIP_VARTYPE_BINARY) );
   SCIP_CALL( SCIPcreateVarBasic(scip, &zvar, "z", -1.0, 1.0, 0.0, SCIP_VARTYPE_CONTINUOUS) );
   SCIP_CALL( SCIPcreateVarBasic(scip, &objvar, "objvar", -SCIPinfinity(scip), SCIPinfinity(scip), 1.0, SCIP_VARTYPE_CONTINUOUS) );

   SCIP_CALL( SCIPaddVar(scip, xvar) );
   SCIP_CALL( SCIPaddVar(scip, yvar) );
   SCIP_CALL( SCIPaddVar(scip, zvar) );
   SCIP_CALL( SCIPaddVar(scip, objvar) );

   /* create quadratic inequality */
   SCIP_CALL( SCIPcreateConsBasicQuadratic(scip, &cons, "quadratic", 0, NULL, NULL, 0, NULL, NULL, NULL, 2.0, 2.0) );

   SCIP_CALL( SCIPaddLinearVarQuadratic(scip, cons, objvar, 1.0) );
   SCIP_CALL( SCIPaddLinearVarQuadratic(scip, cons, zvar, -3.0) );

   SCIP_CALL( SCIPaddSquareCoefQuadratic(scip, cons, xvar, 2.0) );
   SCIP_CALL( SCIPaddSquareCoefQuadratic(scip, cons, yvar, -1.0) );
   SCIP_CALL( SCIPaddSquareCoefQuadratic(scip, cons, zvar, 1.0) );

   SCIP_CALL( SCIPaddBilinTermQuadratic(scip, cons, xvar, yvar, -2.0) );
   SCIP_CALL( SCIPaddBilinTermQuadratic(scip, cons, xvar, zvar, 1.0) );
   SCIP_CALL( SCIPaddBilinTermQuadratic(scip, cons, zvar, yvar, 1.0) );

   SCIP_CALL( SCIPaddCons(scip, cons) );
   SCIP_CALL( SCIPreleaseCons(scip, &cons) );

   /* create knapsack inequality */
   vars[0] = xvar;
   vars[1] = yvar;

   vals[0] = 1;
   vals[1] = 1;

   SCIP_CALL( SCIPcreateConsBasicKnapsack(scip, &cons, "lower", 2, vars, vals, 1) );
   SCIP_CALL( SCIPaddCons(scip, cons) );

   /* check knapsack constraint */
   SCIP_CALL( checkConsKnapsack(scip, cons) );

   SCIP_CALL( SCIPreleaseCons(scip, &cons) );
   SCIP_CALL( SCIPreleaseVar(scip, &xvar) );
   SCIP_CALL( SCIPreleaseVar(scip, &yvar) );
   SCIP_CALL( SCIPreleaseVar(scip, &zvar) );
   SCIP_CALL( SCIPreleaseVar(scip, &objvar) );

   SCIP_CALL( SCIPprintOrigProblem(scip, NULL, NULL, FALSE) );

   return SCIP_OKAY;
}


/** create quadratic problem with one set packing constraint */
static
SCIP_RETCODE createProbSetppc(
   SCIP*                 scip                /**< SCIP instance */
   )
{
   SCIP_CONS* cons;
   SCIP_VAR* xvar;
   SCIP_VAR* yvar;
   SCIP_VAR* zvar;
   SCIP_VAR* objvar;
   SCIP_VAR* vars[3];

   SCIP_CALL( SCIPcreateProbBasic(scip, "qpkktref-setppc") );

   /* create variables */
   SCIP_CALL( SCIPcreateVarBasic(scip, &xvar, "x", 0.0, 1.0, 0.0, SCIP_VARTYPE_BINARY) );
   SCIP_CALL( SCIPcreateVarBasic(scip, &yvar, "y", 0.0, 1.0, 0.0, SCIP_VARTYPE_BINARY) );
   SCIP_CALL( SCIPcreateVarBasic(scip, &zvar, "z", 0.0, 1.0, 0.0, SCIP_VARTYPE_BINARY) );
   SCIP_CALL( SCIPcreateVarBasic(scip, &objvar, "objvar", -SCIPinfinity(scip), SCIPinfinity(scip), 1.0, SCIP_VARTYPE_CONTINUOUS) );

   SCIP_CALL( SCIPaddVar(scip, xvar) );
   SCIP_CALL( SCIPaddVar(scip, yvar) );
   SCIP_CALL( SCIPaddVar(scip, zvar) );
   SCIP_CALL( SCIPaddVar(scip, objvar) );

   /* create quadratic inequality */
   SCIP_CALL( SCIPcreateConsBasicQuadratic(scip, &cons, "quadratic", 0, NULL, NULL, 0, NULL, NULL, NULL, 2.0, 2.0) );

   SCIP_CALL( SCIPaddLinearVarQuadratic(scip, cons, objvar, 1.0) );
   SCIP_CALL( SCIPaddLinearVarQuadratic(scip, cons, zvar, -3.0) );

   SCIP_CALL( SCIPaddSquareCoefQuadratic(scip, cons, xvar, 2.0) );
   SCIP_CALL( SCIPaddSquareCoefQuadratic(scip, cons, yvar, -1.0) );
   SCIP_CALL( SCIPaddSquareCoefQuadratic(scip, cons, zvar, 1.0) );

   SCIP_CALL( SCIPaddBilinTermQuadratic(scip, cons, xvar, yvar, -2.0) );
   SCIP_CALL( SCIPaddBilinTermQuadratic(scip, cons, xvar, zvar, 1.0) );
   SCIP_CALL( SCIPaddBilinTermQuadratic(scip, cons, zvar, yvar, 1.0) );

   SCIP_CALL( SCIPaddCons(scip, cons) );
   SCIP_CALL( SCIPreleaseCons(scip, &cons) );

   /* create set packing inequality */
   vars[0] = xvar;
   vars[1] = yvar;
   vars[2] = zvar;

   SCIP_CALL( SCIPcreateConsBasicSetpack(scip, &cons, "lower", 3, vars) );
   SCIP_CALL( SCIPaddCons(scip, cons) );

   SCIP_CALL( SCIPreleaseCons(scip, &cons) );
   SCIP_CALL( SCIPreleaseVar(scip, &xvar) );
   SCIP_CALL( SCIPreleaseVar(scip, &yvar) );
   SCIP_CALL( SCIPreleaseVar(scip, &zvar) );
   SCIP_CALL( SCIPreleaseVar(scip, &objvar) );

   SCIP_CALL( SCIPprintOrigProblem(scip, NULL, NULL, FALSE) );

   return SCIP_OKAY;
}


/** create quadratic problem with one logicor constraint */
static
SCIP_RETCODE createProbLogicor(
   SCIP*                 scip                /**< SCIP instance */
   )
{
   SCIP_CONS* cons;
   SCIP_VAR* xvar;
   SCIP_VAR* yvar;
   SCIP_VAR* zvar;
   SCIP_VAR* objvar;
   SCIP_VAR* vars[3];

   SCIP_CALL( SCIPcreateProbBasic(scip, "qpkktref-logicor") );

   /* create variables */
   SCIP_CALL( SCIPcreateVarBasic(scip, &xvar, "x", 0.0, 1.0, 0.0, SCIP_VARTYPE_BINARY) );
   SCIP_CALL( SCIPcreateVarBasic(scip, &yvar, "y", 0.0, 1.0, 0.0, SCIP_VARTYPE_BINARY) );
   SCIP_CALL( SCIPcreateVarBasic(scip, &zvar, "z", 0.0, 1.0, 0.0, SCIP_VARTYPE_BINARY) );
   SCIP_CALL( SCIPcreateVarBasic(scip, &objvar, "objvar", -SCIPinfinity(scip), SCIPinfinity(scip), 1.0, SCIP_VARTYPE_CONTINUOUS) );

   SCIP_CALL( SCIPaddVar(scip, xvar) );
   SCIP_CALL( SCIPaddVar(scip, yvar) );
   SCIP_CALL( SCIPaddVar(scip, zvar) );
   SCIP_CALL( SCIPaddVar(scip, objvar) );

   /* create quadratic inequality */
   SCIP_CALL( SCIPcreateConsBasicQuadratic(scip, &cons, "quadratic", 0, NULL, NULL, 0, NULL, NULL, NULL, 2.0, 2.0) );

   SCIP_CALL( SCIPaddLinearVarQuadratic(scip, cons, objvar, -1.0) );
   SCIP_CALL( SCIPaddLinearVarQuadratic(scip, cons, zvar, -1.0) );

   SCIP_CALL( SCIPaddSquareCoefQuadratic(scip, cons, xvar, 1.0) );
   SCIP_CALL( SCIPaddSquareCoefQuadratic(scip, cons, yvar, -1.0) );
   SCIP_CALL( SCIPaddSquareCoefQuadratic(scip, cons, zvar, -2.0) );

   SCIP_CALL( SCIPaddBilinTermQuadratic(scip, cons, xvar, yvar, 2.0) );
   SCIP_CALL( SCIPaddBilinTermQuadratic(scip, cons, xvar, zvar, 1.0) );
   SCIP_CALL( SCIPaddBilinTermQuadratic(scip, cons, zvar, yvar, 1.0) );

   SCIP_CALL( SCIPaddCons(scip, cons) );
   SCIP_CALL( SCIPreleaseCons(scip, &cons) );

   /* create logicor inequality */
   vars[0] = xvar;
   vars[1] = yvar;
   vars[2] = zvar;

   SCIP_CALL( SCIPcreateConsBasicLogicor(scip, &cons, "lower", 3, vars) );
   SCIP_CALL( SCIPaddCons(scip, cons) );

   SCIP_CALL( SCIPreleaseCons(scip, &cons) );
   SCIP_CALL( SCIPreleaseVar(scip, &xvar) );
   SCIP_CALL( SCIPreleaseVar(scip, &yvar) );
   SCIP_CALL( SCIPreleaseVar(scip, &zvar) );
   SCIP_CALL( SCIPreleaseVar(scip, &objvar) );

   SCIP_CALL( SCIPprintOrigProblem(scip, NULL, NULL, FALSE) );

   return SCIP_OKAY;
}


/** create quadratic problem with one varbound */
static
SCIP_RETCODE createProbVarbound(
   SCIP*                 scip                /**< SCIP instance */
   )
{
   SCIP_CONS* cons;
   SCIP_VAR* xvar;
   SCIP_VAR* yvar;
   SCIP_VAR* zvar;
   SCIP_VAR* objvar;

   SCIP_CALL( SCIPcreateProbBasic(scip, "qpkktref-varbound") );

   SCIP_CALL( SCIPsetObjsense(scip, SCIP_OBJSENSE_MAXIMIZE) );

   /* create variables */
   SCIP_CALL( SCIPcreateVarBasic(scip, &xvar, "x", 0.0, 1.0, 0.0, SCIP_VARTYPE_BINARY) );
   SCIP_CALL( SCIPcreateVarBasic(scip, &yvar, "y", 0.0, 1.0, 0.0, SCIP_VARTYPE_BINARY) );
   SCIP_CALL( SCIPcreateVarBasic(scip, &zvar, "z", -1.0, 1.0, 0.0, SCIP_VARTYPE_CONTINUOUS) );
   SCIP_CALL( SCIPcreateVarBasic(scip, &objvar, "objvar", -SCIPinfinity(scip), SCIPinfinity(scip), 1.0, SCIP_VARTYPE_CONTINUOUS) );

   SCIP_CALL( SCIPaddVar(scip, xvar) );
   SCIP_CALL( SCIPaddVar(scip, yvar) );
   SCIP_CALL( SCIPaddVar(scip, zvar) );
   SCIP_CALL( SCIPaddVar(scip, objvar) );

   /* create quadratic inequality */
   SCIP_CALL( SCIPcreateConsBasicQuadratic(scip, &cons, "quadratic", 0, NULL, NULL, 0, NULL, NULL, NULL, 2.0, 2.0) );

   SCIP_CALL( SCIPaddLinearVarQuadratic(scip, cons, objvar, 1.0) );
   SCIP_CALL( SCIPaddLinearVarQuadratic(scip, cons, zvar, -3.0) );

   SCIP_CALL( SCIPaddSquareCoefQuadratic(scip, cons, xvar, 2.0) );
   SCIP_CALL( SCIPaddSquareCoefQuadratic(scip, cons, yvar, -1.0) );
   SCIP_CALL( SCIPaddSquareCoefQuadratic(scip, cons, zvar, 1.0) );

   SCIP_CALL( SCIPaddBilinTermQuadratic(scip, cons, xvar, yvar, -2.0) );
   SCIP_CALL( SCIPaddBilinTermQuadratic(scip, cons, xvar, zvar, 1.0) );
   SCIP_CALL( SCIPaddBilinTermQuadratic(scip, cons, zvar, yvar, 1.0) );

   SCIP_CALL( SCIPaddCons(scip, cons) );
   SCIP_CALL( SCIPreleaseCons(scip, &cons) );

   /* create varbound inequality */
   SCIP_CALL( SCIPcreateConsBasicVarbound(scip, &cons, "lower", zvar, xvar, 1.0, -1.0, 1.0) );
   SCIP_CALL( SCIPaddCons(scip, cons) );

   SCIP_CALL( SCIPreleaseCons(scip, &cons) );
   SCIP_CALL( SCIPreleaseVar(scip, &xvar) );
   SCIP_CALL( SCIPreleaseVar(scip, &yvar) );
   SCIP_CALL( SCIPreleaseVar(scip, &zvar) );
   SCIP_CALL( SCIPreleaseVar(scip, &objvar) );

   SCIP_CALL( SCIPprintOrigProblem(scip, NULL, NULL, FALSE) );

   return SCIP_OKAY;
}


/** create quadratic problem */
static
SCIP_RETCODE createProb(
   SCIP*                 scip,               /**< SCIP instance */
   int                   instance            /**< instance to create */
   )
{
   switch( instance )
   {
      case 0:
         SCIP_CALL( createProbLinear1(scip) );
         break;
      case 1:
         SCIP_CALL( createProbLinear2(scip) );
         break;
      case 2:
         SCIP_CALL( createProbKnapsack(scip) );
         break;
      case 3:
         SCIP_CALL( createProbSetppc(scip) );
        break;
      case 4:
         SCIP_CALL( createProbLogicor(scip) );
        break;
      case 5:
         SCIP_CALL( createProbVarbound(scip) );
        break;
      default:
         SCIPerrorMessage("unknown instance number\n");
         return SCIP_INVALIDDATA;
   }

   return SCIP_OKAY;
}

/** include QP settings (KKT conditions are not added, all presolvers are turned off) */
static
SCIP_RETCODE includeQPSettings(
   SCIP*                 scip                /**< SCIP instance */
   )
{
   assert( scip != NULL );

   /* include default plugins */
   SCIP_CALL( SCIPincludeDefaultPlugins(scip) );

   SCIP_CALL( SCIPsetIntParam(scip, "presolving/qpkktref/maxrounds", 0) );

   return SCIP_OKAY;
}


/** include KKT settings (KKT conditions are added, all presolvers are turned off) */
static
SCIP_RETCODE includeKKTSettings(
   SCIP*                 scip                /**< SCIP instance */
   )
{
   assert( scip != NULL );

   SCIP_CALL( includeQPSettings(scip) );

   SCIP_CALL( SCIPsetIntParam(scip, "presolving/qpkktref/maxrounds", -1) );

   /* allow variables to be unbounded */
   SCIP_CALL( SCIPsetBoolParam(scip, "presolving/qpkktref/updatequadbounded", FALSE) );

   /* allow binary variables */
   SCIP_CALL( SCIPsetBoolParam(scip, "presolving/qpkktref/addkktbinary", TRUE) );

   SCIP_CALL( SCIPsetIntParam(scip, "constraints/linear/maxprerounds", 0) );
   SCIP_CALL( SCIPsetIntParam(scip, "constraints/knapsack/maxprerounds", 0) );
   SCIP_CALL( SCIPsetIntParam(scip, "constraints/varbound/maxprerounds", 0) );
   SCIP_CALL( SCIPsetIntParam(scip, "constraints/logicor/maxprerounds", 0) );
   SCIP_CALL( SCIPsetIntParam(scip, "constraints/setppc/maxprerounds", 0) );
   SCIP_CALL( SCIPsetIntParam(scip, "constraints/quadratic/maxprerounds", 0) );

   return SCIP_OKAY;
}


/** run unittest */
static
SCIP_RETCODE runUnittest(void)
{
   SCIP* scip1 = NULL;
   SCIP* scip2 = NULL;
   int j;

   /* output stuff for automatic unittest evaluation */
   printf("@01 unittest-qpkktref ===========\n");
   printf("=opt=  unittest-qpkktref 0\n\n");

   for (j = 0; j < 6; ++j)
   {
      SCIP_Bool equal = FALSE;

      /* initialize SCIP */
      SCIP_CALL( SCIPcreate(&scip1) );
      SCIP_CALL( SCIPcreate(&scip2) );

      /* include QP settings to scip1 (KKT conditions are not added, all presolvers are turned off) */
      SCIP_CALL( includeQPSettings(scip1) );

      /* include KKT settings to scip2 (KKT conditions are added, all presolvers are turned off) */
      SCIP_CALL( includeKKTSettings(scip2) );

      /* create problem */
      SCIP_CALL( createProb(scip1, j) );
      SCIP_CALL( createProb(scip2, j) );

      /* solve */
      SCIP_CALL( SCIPsolve(scip1) );
      SCIP_CALL( SCIPsolve(scip2) );
      /*
      SCIP_CALL( SCIPprintBestSol(scip1, NULL, FALSE) );
      SCIP_CALL( SCIPprintStatistics(scip1, NULL) );
      */
      /*SCIP_CALL( SCIPwriteTrafoProblem(scip2, "trafounittestQP.lp", NULL, FALSE ) );*/

      if ( SCIPisFeasEQ(scip1, SCIPgetPrimalbound(scip1), SCIPgetPrimalbound(scip2) ) )
         equal = TRUE;

      /* free transformed problem */
      SCIP_CALL( SCIPfreeTransform(scip1) );
      SCIP_CALL( SCIPfreeTransform(scip2) );

      /* free SCIP */
      SCIP_CALL( SCIPfree(&scip1) );
      SCIP_CALL( SCIPfree(&scip2) );

      /* check for memory leaks */
      BMScheckEmptyMemory();

      if ( ! equal )
      {
         SCIPerrorMessage("Optimal solution of original problem is not equal to optimal solution of reformulated problem.\n");
         return SCIP_ERROR;
      }
   }

   printf("Unit test for KKT-reformulation passed\n");

   /* for automatic testing output the following */
   printf("SCIP Status        : all tests passed\n");
   printf("Ignore the following:\n");
   printf("  solving          : 0.00\n");
   printf("  nodes (total)    : 0\n");
   printf("  Primal Bound     : 0.0\n");
   printf("  Dual Bound       : 0.0\n");

   return SCIP_OKAY;
}

/** main function */
int
main(
   int                        argc,
   char**                     argv
   )
{
   SCIP_RETCODE retcode;

   retcode = runUnittest();

   if ( retcode != SCIP_OKAY )
   {
      SCIPprintError(retcode);
      return -1;
   }

   return 0;
}
